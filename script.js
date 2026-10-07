// =====================================================================
// Pendule inverse sur chariot - simulation temps reel
// Convention : x vers la droite, theta = 0 pendule vertical VERS LE HAUT,
// theta > 0 quand le pendule penche vers la droite.
// =====================================================================

const clamp = (v, lo, hi) => Math.max(lo, Math.min(hi, v));

// Les textes numeriques ne sont rafraichis que 4 fois par seconde :
// lisibles, sans chiffres qui scintillent a chaque image.
const TEXT_REFRESH_MS = 250;
const textCache = new Map();
const steady = (key, make) => {
    const now = Date.now();
    const c = textCache.get(key);
    if (c && now - c.t < TEXT_REFRESH_MS) return c.s;
    const s = make();
    textCache.set(key, { t: now, s });
    return s;
};

const wrapAngle = (a) => {
    const twoPi = 2 * Math.PI;
    a = (a + Math.PI) % twoPi;
    if (a < 0) a += twoPi;
    return a - Math.PI;
};

// ===== MODULE PHYSICS =====
class Physics {
    constructor(params = {}) {
        this.mc = 1.0;      // masse chariot (kg)
        this.mp = 0.3;      // masse pendule (kg)
        this.l = 0.5;       // distance pivot -> masse (m)
        this.b = 0.05;      // frottement chariot (N.s/m)
        this.c = 0.01;      // frottement pivot (N.m.s)
        this.g = 9.81;      // gravite (m/s2)
        this.fMax = 50;     // saturation de l'actionneur (N)
        this.railLimit = 2.2; // demi-longueur du rail (m)
        this.coulomb = 0;   // frottement sec du chariot (N)
        Object.assign(this, params);

        this.state = [0, 0, 0.1, 0]; // [x, x_dot, theta, theta_dot]
        this.t = 0;
        this.dt = 0.005;
        this.lastForce = 0;
    }

    updateParams(params) {
        Object.assign(this, params);
    }

    // Le pendule est une tige homogene de longueur 2l (centre de masse a l du pivot),
    // inertie propre I = mp l^2 / 3, donc longueur equivalente Le = l + I/(mp l) = 4l/3.
    get lEq() { return 4 * this.l / 3; }

    // Equations de Lagrange :
    //  (mc+mp) x'' + mp l cos(th) th'' = F - b x' - Fs sgn(x') + mp l th'^2 sin(th)
    //  cos(th) x'' + Le th''           = g sin(th) - c/(mp l) th'
    derivatives(state, F) {
        const [, x_dot, theta, theta_dot] = state;
        const S = Math.sin(theta);
        const C = Math.cos(theta);

        const a11 = this.mc + this.mp;
        const a12 = this.mp * this.l * C;
        const a21 = C;
        const a22 = this.lEq;

        // Frottement sec lisse (tanh) pour rester integrable
        const dry = this.coulomb * Math.tanh(x_dot / 0.01);
        const v1 = F - this.b * x_dot - dry + this.mp * this.l * theta_dot * theta_dot * S;
        const v2 = this.g * S - (this.c / (this.mp * this.l)) * theta_dot;

        // det = (mc+mp) Le - mp l cos^2 > 0 toujours
        const det = a11 * a22 - a12 * a21;
        const x_ddot = (a22 * v1 - a12 * v2) / det;
        const theta_ddot = (-a21 * v1 + a11 * v2) / det;

        return [x_dot, x_ddot, theta_dot, theta_ddot];
    }

    rk4Step(F) {
        const s = this.state, h = this.dt;
        const k1 = this.derivatives(s, F);
        const k2 = this.derivatives(s.map((v, i) => v + h * k1[i] / 2), F);
        const k3 = this.derivatives(s.map((v, i) => v + h * k2[i] / 2), F);
        const k4 = this.derivatives(s.map((v, i) => v + h * k3[i]), F);
        for (let i = 0; i < 4; i++) {
            s[i] += h * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6;
        }
    }

    eulerStep(F) {
        const d = this.derivatives(this.state, F);
        for (let i = 0; i < 4; i++) this.state[i] += this.dt * d[i];
    }

    // force : moteur (sature a fMax) ; hand : force exterieure de la main (non saturee)
    step(force, integrator = 'rk4', hand = 0) {
        this.lastForce = clamp(force, -this.fMax, this.fMax);
        this.lastHand = hand;
        const F = this.lastForce + hand;

        if (integrator === 'rk4') this.rk4Step(F);
        else this.eulerStep(F);

        this.state[2] = wrapAngle(this.state[2]);

        // Butees du rail : choc inelastique
        if (Math.abs(this.state[0]) > this.railLimit) {
            this.state[0] = Math.sign(this.state[0]) * this.railLimit;
            this.state[1] = 0;
        }

        this.t += this.dt;
    }

    reset(initialState = [0, 0, 0.1, 0]) {
        this.state = [...initialState];
        this.t = 0;
        this.lastForce = 0;
    }

    // Impulsion horizontale (N.s) appliquee au chariot
    applyImpulse(impulse) {
        this.state[1] += impulse / this.mc;
    }

    // Bilan des forces a l'instant courant (N), repere x vers la droite, y vers le haut.
    // La force du pivot sur la tige vient de Newton applique au centre de masse G :
    //   R = mp aG - mp g   avec   G = (x + l sin th, l cos th)
    forces(F) {
        const [, x_dot, theta, theta_dot] = this.state;
        const hand = this.lastHand || 0;
        const [, x_ddot, , theta_ddot] = this.derivatives(this.state, F + hand);
        const S = Math.sin(theta), C = Math.cos(theta);
        const aGx = x_ddot + this.l * (C * theta_ddot - S * theta_dot ** 2);
        const aGy = -this.l * (S * theta_ddot + C * theta_dot ** 2);
        const Rx = this.mp * aGx;
        const Ry = this.mp * aGy + this.mp * this.g;
        return {
            motor: F,
            hand,
            weightPole: this.mp * this.g,
            weightCart: this.mc * this.g,
            pivot: [Rx, Ry],                       // force du chariot sur la tige
            normal: this.mc * this.g + Ry,         // reaction du rail sur le chariot
            friction: -this.b * x_dot - this.coulomb * Math.tanh(x_dot / 0.01)
        };
    }

    setState(x, x_dot, theta, theta_dot) {
        this.state = [x, x_dot, theta, theta_dot];
    }

    get x() { return this.state[0]; }
    get x_dot() { return this.state[1]; }
    get theta() { return this.state[2]; }
    get theta_dot() { return this.state[3]; }
}

// ===== MODULE CONTROLLER =====
// Retour d'etat avec integrale sur l'angle :
//   F = Kp th + Ki int(th) + Kd th' + Kx (x - x_ref) + Kv x'
// Tous les gains sont positifs pour un pendule inverse (il faut deplacer
// le chariot SOUS la masse, donc dans le sens de l'inclinaison).
class PIDController {
    constructor(gains = {}) {
        this.kp = 82.6; this.ki = 0; this.kd = 17.4;
        this.kx = 12.6; this.kv = 15.7;
        Object.assign(this, gains);
        this.integral = 0;
        this.maxIntegral = 0.5;
        this.catchAngle = 0.4;   // |theta| sous lequel on passe en stabilisation (23 deg)
        this.dropAngle = 0.8;    // |theta| au-dela duquel on repasse en redressement
        this.swingUp = true;     // redressement automatique depuis la position basse
        this.swingGain = 8;      // gain du pompage d'energie
        this.swingAccel = 1.2;   // acceleration max du chariot en redressement (en g)
        this.swingKx = 1.5;      // rappel du chariot vers le centre pendant le redressement
        this.swingKv = 3.0;
        this.swingMargin = 0.2;  // vise un peu plus que l'energie du sommet pour l'atteindre
        this.phase = 'swing';    // 'swing' (redressement) | 'balance' (stabilisation)
        this.plant = null;       // parametres physiques (masses, l, g) pour le redressement
    }

    setGains(gains) {
        Object.assign(this, gains);
    }

    compute(state, xRef, dt) {
        const [x, x_dot, theta, theta_dot] = state;
        const th = wrapAngle(theta);

        // Machine a etats avec hysteresis : redressement <-> stabilisation
        if (this.phase === 'balance' && Math.abs(th) > this.dropAngle) {
            this.phase = 'swing';
            this.integral = 0;
        } else if (this.phase === 'swing' && Math.abs(th) < this.catchAngle) {
            this.phase = 'balance';
            this.integral = 0;
        }

        if (this.phase === 'swing') {
            return this.swingUp && this.plant ? this.swingForce(state, xRef) : 0;
        }

        this.integral = clamp(this.integral + th * dt, -this.maxIntegral, this.maxIntegral);

        return this.kp * th + this.ki * this.integral + this.kd * theta_dot
             + this.kx * (x - xRef) + this.kv * x_dot;
    }

    // Redressement par pompage d'energie (Astrom & Furuta).
    // E = 1/2 J th'^2 + mp g l (cos th - 1) vaut 0 en haut et -2 mp g l en bas.
    // dE/dt = -mp l cos(th) th' x'' : on accelere le chariot pour faire monter E vers 0.
    swingForce(state, xRef) {
        const [x, x_dot, theta, theta_dot] = state;
        const { mc, mp, l, g } = this.plant;
        const J = mp * l * (4 * l / 3);                 // inertie de la tige autour du pivot
        const E = 0.5 * J * theta_dot ** 2 + mp * g * l * (Math.cos(theta) - 1);
        const E0 = 2 * mp * g * l;

        // Au repos parfait en bas, sign(0) = 0 : on donne un premier coup
        const dir = Math.abs(theta_dot) > 1e-3 ? Math.sign(theta_dot * Math.cos(theta)) : 1;
        const aMax = this.swingAccel * g;
        let a = clamp(this.swingGain * g * (E / E0 - this.swingMargin) * dir, -aMax, aMax);

        // Rappel doux du chariot vers le centre pour rester sur le rail
        a += -this.swingKx * (x - xRef) - this.swingKv * x_dot;

        return (mc + mp) * a;
    }

    reset() {
        this.integral = 0;
        this.phase = 'swing';
    }

    // Auto-reglage par placement de poles (Ackermann) sur le modele linearise.
    // speed > 1 = reponse plus rapide (et forces plus grandes).
    static autoTune(p, speed = 1) {
        const { mc, mp, l, b, c, g } = p;
        const M = mc + mp, Le = 4 * l / 3;
        const cp = c / (mp * l);
        const det = M * Le - mp * l;
        const A = [
            [0, 1, 0, 0],
            [0, -b * Le / det, -mp * l * g / det, mp * l * cp / det],
            [0, 0, 0, 1],
            [0, b / det, M * g / det, -M * cp / det]
        ];
        const B = [0, Le / det, 0, -1 / det];
        const mul = (X, Y) => X.map(r => Y[0].map((_, j) => r.reduce((s, v, k) => s + v * Y[k][j], 0)));
        const mv = (X, v) => X.map(r => r.reduce((s, a, k) => s + a * v[k], 0));
        const I = [0, 1, 2, 3].map(i => [0, 1, 2, 3].map(j => (i === j ? 1 : 0)));

        // 2 poles lents (position) + 2 poles rapides (angle) lies a la pulsation propre
        const w0 = Math.sqrt(Math.abs(g) * M / det) || 1; // pole instable en boucle ouverte
        const poles = [-1.8 * speed, -2.4 * speed, -1.4 * w0 * speed, -1.6 * w0 * speed];

        // Polynome caracteristique desire
        let coef = [1];
        for (const pk of poles) {
            const n = [...coef, 0];
            for (let i = 1; i < n.length; i++) n[i] -= pk * coef[i - 1];
            coef = n;
        }

        // phi(A) = A^4 + a1 A^3 + a2 A^2 + a3 A + a4 I
        let phi = I.map(r => r.map(v => v * coef[4]));
        let Ak = I;
        for (let k = 1; k <= 4; k++) {
            Ak = mul(Ak, A);
            phi = phi.map((r, i) => r.map((v, j) => v + coef[4 - k] * Ak[i][j]));
        }

        // Commandabilite : on resout C' y = e4 (Gauss-Jordan), puis K = y' phi(A)
        const cols = [B];
        for (let k = 1; k < 4; k++) cols.push(mv(A, cols[k - 1]));
        const G = [0, 1, 2, 3].map(i => [...cols[i], i === 3 ? 1 : 0]);
        for (let i = 0; i < 4; i++) {
            let piv = i;
            for (let r = i + 1; r < 4; r++) if (Math.abs(G[r][i]) > Math.abs(G[piv][i])) piv = r;
            [G[i], G[piv]] = [G[piv], G[i]];
            for (let r = 0; r < 4; r++) {
                if (r === i) continue;
                const f = G[r][i] / G[i][i];
                for (let k = i; k < 5; k++) G[r][k] -= f * G[i][k];
            }
        }
        const y = G.map((r, i) => r[4] / r[i]);
        const K = [0, 1, 2, 3].map(j => y.reduce((s, v, k) => s + v * phi[k][j], 0));

        // F = -K s
        return { kx: -K[0], kv: -K[1], kp: -K[2], ki: 0, kd: -K[3] };
    }

    static presetSpeed(type) {
        return { soft: 0.7, nominal: 1.0, aggressive: 1.4 }[type] || 1.0;
    }
}

// ===== MODULE CAPTEURS =====
// Un vrai banc ne mesure que x et theta (encodeurs), avec bruit et quantification.
// Les vitesses sont estimees par difference finie filtree (passe-bas 1er ordre).
class Sensors {
    constructor() {
        this.enabled = true;
        this.noiseTheta = 0.002;              // ecart-type du bruit angulaire (rad)
        this.qTheta = 2 * Math.PI / 4096;     // encodeur 4096 points/tour
        this.qX = 0.0001;                     // encodeur lineaire 0.1 mm
        this.cutoffHz = 25;                   // filtre des vitesses estimees
        this.reset();
    }

    reset() {
        this.prev = null;
        this.vx = 0;
        this.vth = 0;
    }

    static gauss() {
        let u = 0, v = 0;
        while (u === 0) u = Math.random();
        while (v === 0) v = Math.random();
        return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
    }

    read(state, Ts) {
        if (!this.enabled) return [...state];

        const [x, , th] = state;
        const q = (v, step) => Math.round(v / step) * step;
        const xm = q(x + Sensors.gauss() * this.noiseTheta * 0.5, this.qX);
        const thm = q(th + Sensors.gauss() * this.noiseTheta, this.qTheta);

        if (this.prev) {
            const a = Ts / (Ts + 1 / (2 * Math.PI * this.cutoffHz));
            this.vx += a * ((xm - this.prev[0]) / Ts - this.vx);
            this.vth += a * (wrapAngle(thm - this.prev[1]) / Ts - this.vth);
        }
        this.prev = [xm, thm];
        return [xm, this.vx, thm, this.vth];
    }
}

// ===== MODULE UI =====
const SLIDERS = ['kp', 'ki', 'kd', 'kx', 'kv', 'mc', 'mp', 'l', 'b', 'c', 'g', 'fmax', 'dt',
                 'noise', 'tau', 'ts', 'fs'];

class UI {
    constructor() {
        this.mode = 'strict'; // 'strict' (stabilisation) ou 'project' (suivi de position)
        this.$ = (id) => document.getElementById(id);

        this.playBtn = this.$('play-btn');
        this.pauseBtn = this.$('pause-btn');
        this.resetBtn = this.$('reset-btn');
        this.impulseBtn = this.$('impulse-btn');
        this.strictModeBtn = this.$('strict-mode');
        this.projectModeBtn = this.$('project-mode');
        this.pidPreset = this.$('pid-preset');
        this.autoTuneBtn = this.$('auto-tune');
        this.integratorSelect = this.$('integrator');
        this.scenarioSelect = this.$('scenario');

        SLIDERS.forEach(name => {
            const slider = this.$(`${name}-slider`);
            slider.addEventListener('input', () => {
                this.refreshValue(name);
                if (['kp', 'ki', 'kd', 'kx', 'kv'].includes(name)) this.pidPreset.value = 'custom';
            });
            this.refreshValue(name);
        });

        this.playBtn.addEventListener('click', () => this.onPlay());
        this.pauseBtn.addEventListener('click', () => this.onPause());
        this.resetBtn.addEventListener('click', () => this.onReset());
        this.impulseBtn.addEventListener('click', () => this.onImpulse());
        this.strictModeBtn.addEventListener('click', () => { this.setMode('strict'); this.onModeChange(); });
        this.projectModeBtn.addEventListener('click', () => { this.setMode('project'); this.onModeChange(); });
        this.autoTuneBtn.addEventListener('click', () => this.onAutoTune());
        this.pidPreset.addEventListener('change', () => this.onPresetChange());
        this.scenarioSelect.addEventListener('change', () => this.onScenarioChange());
    }

    refreshValue(name) {
        const v = parseFloat(this.$(`${name}-slider`).value);
        const unit = { dt: ' ms', ts: ' ms', tau: ' ms', noise: ' mrad', fs: ' N' }[name];
        this.$(`${name}-value`).textContent = unit ? `${+v.toFixed(2)}${unit}` : v.toFixed(2);
    }

    getRealismParams() {
        const v = (n) => parseFloat(this.$(`${n}-slider`).value);
        return {
            enabled: this.$('realism').checked,
            noise: v('noise') / 1000,   // mrad -> rad
            tau: v('tau') / 1000,       // ms -> s
            ts: v('ts') / 1000,         // ms -> s
            fs: v('fs')
        };
    }

    setMode(mode) {
        this.mode = mode;
        this.strictModeBtn.classList.toggle('active', mode === 'strict');
        this.projectModeBtn.classList.toggle('active', mode === 'project');
    }

    // Callbacks branches par l'application
    onPlay() {}
    onPause() {}
    onReset() {}
    onImpulse() {}
    onAutoTune() {}
    onPresetChange() {}
    onScenarioChange() {}
    onModeChange() {}

    updateIndicators(m) {
        // Rafraichi 4 fois par seconde seulement (evite chiffres et couleurs qui clignotent)
        const now = Date.now();
        if (this.lastIndicators && now - this.lastIndicators < TEXT_REFRESH_MS) return;
        this.lastIndicators = now;

        this.$('stability-time').textContent = m.stabilityTime !== null ? `${m.stabilityTime.toFixed(2)} s` : '-- s';
        this.$('max-angle').textContent = `${m.maxAngle.toFixed(3)} rad`;
        this.$('rise-time').textContent = m.riseTime !== null ? `${m.riseTime.toFixed(2)} s` : '-- s';
        this.$('steady-error').textContent = m.steadyError !== null ? `${m.steadyError.toFixed(1)} %` : '-- %';
        this.$('system-status').textContent = m.systemStatus;
        this.$('force-value').textContent = `${m.force.toFixed(1)} N`;

        this.setStatus('stability', m.stabilityStatus);
        this.setStatus('angle', m.angleStatus);
        this.setStatus('rise', m.riseStatus);
        this.setStatus('error', m.errorStatus);
        this.setStatus('status', m.overallStatus);

        const projectOnly = this.mode !== 'project';
        this.$('rise-indicator').style.opacity = projectOnly ? 0.4 : 1;
        this.$('error-indicator').style.opacity = projectOnly ? 0.4 : 1;
    }

    setStatus(type, status) {
        const el = this.$(`${type}-indicator`);
        const cls = 'indicator' + (status ? ` ${status}` : '');
        if (el.className !== cls) el.className = cls;
    }

    getGains() {
        const v = (n) => parseFloat(this.$(`${n}-slider`).value);
        return { kp: v('kp'), ki: v('ki'), kd: v('kd'), kx: v('kx'), kv: v('kv') };
    }

    setGains(gains) {
        for (const k of ['kp', 'ki', 'kd', 'kx', 'kv']) {
            const slider = this.$(`${k}-slider`);
            // Elargit la plage si le gain calcule depasse le max du slider
            if (gains[k] > parseFloat(slider.max)) slider.max = Math.ceil(gains[k] * 1.5);
            slider.value = gains[k];
            this.refreshValue(k);
        }
    }

    getSystemParams() {
        const v = (n) => parseFloat(this.$(`${n}-slider`).value);
        return { mc: v('mc'), mp: v('mp'), l: v('l'), b: v('b'), c: v('c'), g: v('g'), fMax: v('fmax') };
    }

    getSimulationParams() {
        return {
            dt: parseFloat(this.$('dt-slider').value) / 1000, // ms -> s
            integrator: this.integratorSelect.value
        };
    }
}

// ===== MODULE ANIMATION AVEC INTERACTIVITE =====
const FORCE_SCALE = 4;     // px par newton
const FORCE_MAX_PX = 220;  // au-dela, la fleche est tronquee (la valeur reste affichee)
const FORCE_COLORS = {
    weight: '#d35400',   // poids
    motor: '#27ae60',    // moteur
    hand: '#2c3e50',     // main (glisser le chariot)
    normal: '#2980b9',   // reaction du rail
    pivot: '#16a085',    // force du pivot
    friction: '#7f8c8d'  // frottements
};

class Animation {
    constructor(canvasId) {
        this.canvas = document.getElementById(canvasId);
        this.ctx = this.canvas.getContext('2d');

        this.scale = 150; // pixels par metre
        this.centerX = this.canvas.width / 2;
        this.cartY = this.canvas.height * 0.5; // place pour la tige en haut comme en bas

        this.cartWidth = 70;
        this.cartHeight = 30;
        this.wheelRadius = 8;
        this.railY = this.cartY + this.cartHeight / 2 + this.wheelRadius;

        this.isDragging = false;
        this.dragTarget = null; // 'cart' | 'pendulum'
        this.mouseX = -1;
        this.mouseY = -1;

        this.onCartDrag = null;
        this.onPendulumDrag = null;
        this.onDragEnd = null;

        this.currentX = 0;
        this.currentTheta = 0;
        this.currentL = 0.5;

        this.setupInteraction();
    }

    setupInteraction() {
        const c = this.canvas;
        c.addEventListener('mousedown', (e) => this.start(this.pos(e)));
        c.addEventListener('mousemove', (e) => this.move(this.pos(e)));
        c.addEventListener('mouseup', () => this.end());
        c.addEventListener('mouseleave', () => this.end());
        c.addEventListener('touchstart', (e) => { e.preventDefault(); this.start(this.pos(e.touches[0])); }, { passive: false });
        c.addEventListener('touchmove', (e) => { e.preventDefault(); this.move(this.pos(e.touches[0])); }, { passive: false });
        c.addEventListener('touchend', (e) => { e.preventDefault(); this.end(); }, { passive: false });
    }

    // Coordonnees souris -> coordonnees canvas (tient compte du redimensionnement CSS)
    pos(e) {
        const r = this.canvas.getBoundingClientRect();
        return {
            x: (e.clientX - r.left) * this.canvas.width / r.width,
            y: (e.clientY - r.top) * this.canvas.height / r.height
        };
    }

    cartPx() { return this.centerX + this.currentX * this.scale; }

    // Extremite de la tige (longueur totale 2l)
    bobPx() {
        const L = 2 * this.currentL * this.scale;
        return {
            x: this.cartPx() + L * Math.sin(this.currentTheta),
            y: this.cartY - L * Math.cos(this.currentTheta)
        };
    }

    hitCart(x, y) {
        const cx = this.cartPx();
        return Math.abs(x - cx) <= this.cartWidth / 2 && Math.abs(y - this.cartY) <= this.cartHeight / 2;
    }

    hitBob(x, y) {
        const b = this.bobPx();
        return (x - b.x) ** 2 + (y - b.y) ** 2 <= 18 * 18;
    }

    start({ x, y }) {
        if (this.hitBob(x, y)) this.dragTarget = 'pendulum';
        else if (this.hitCart(x, y)) this.dragTarget = 'cart';
        else return;
        this.isDragging = true;
        this.canvas.style.cursor = 'grabbing';
        this.mouseX = x; this.mouseY = y;
    }

    move({ x, y }) {
        if (!this.isDragging) {
            this.canvas.style.cursor = (this.hitBob(x, y) || this.hitCart(x, y)) ? 'grab' : 'default';
            this.mouseX = x; this.mouseY = y;
            return;
        }
        if (this.dragTarget === 'cart' && this.onCartDrag) {
            this.onCartDrag((x - this.mouseX) / this.scale);
        } else if (this.dragTarget === 'pendulum' && this.onPendulumDrag) {
            this.onPendulumDrag(Math.atan2(x - this.cartPx(), this.cartY - y));
        }
        this.mouseX = x; this.mouseY = y;
    }

    end() {
        if (this.isDragging && this.onDragEnd) this.onDragEnd();
        this.isDragging = false;
        this.dragTarget = null;
        this.canvas.style.cursor = 'default';
    }

    draw(x, theta, l, force, targetX = null, forces = null) {
        this.currentX = x;
        this.currentTheta = theta;
        this.currentL = l;
        // Echelle reduite si la tige ne tient pas au-dessus du chariot
        this.scale = Math.min(150, (this.cartY - 30) / (2 * l));

        const ctx = this.ctx;
        ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);

        this.drawRail();
        if (targetX !== null) this.drawTarget(this.centerX + targetX * this.scale);

        const cartX = this.cartPx();
        if (!forces) this.drawForce(cartX, force);
        this.drawCart(cartX);
        this.drawAngle(cartX, theta);
        this.drawPole(cartX);
        if (forces) {
            this.drawForceVectors(cartX, theta, forces);
            this.drawForceLegend();
        }
        this.drawValues(x, theta, force);
        this.drawInstructions();
    }

    // Longueur affichee proportionnelle a la racine de |F| : petites et grandes forces restent lisibles
    // Longueur proportionnelle a la force : FORCE_SCALE px par newton
    static arrowLength(F) {
        return Math.min(FORCE_MAX_PX, FORCE_SCALE * Math.abs(F));
    }

    // Fleche partant de (x0, y0) dans la direction (ux, uy) (repere canvas), de longueur len
    arrow(x0, y0, ux, uy, len, color, label, labelSide = 1) {
        if (len < 1.5) return;
        const ctx = this.ctx;
        const x1 = x0 + ux * len, y1 = y0 + uy * len;
        ctx.strokeStyle = ctx.fillStyle = color;
        ctx.lineWidth = 2.5;
        ctx.beginPath();
        ctx.moveTo(x0, y0);
        ctx.lineTo(x1, y1);
        ctx.stroke();
        const nx = -uy, ny = ux;
        ctx.beginPath();
        ctx.moveTo(x1, y1);
        ctx.lineTo(x1 - 9 * ux + 4.5 * nx, y1 - 9 * uy + 4.5 * ny);
        ctx.lineTo(x1 - 9 * ux - 4.5 * nx, y1 - 9 * uy - 4.5 * ny);
        ctx.closePath();
        ctx.fill();
        if (label) {
            ctx.font = 'bold 12px sans-serif';
            ctx.textAlign = 'center';
            ctx.textBaseline = 'middle';
            const lx = x1 + ux * 14 + nx * 12 * labelSide, ly = y1 + uy * 14 + ny * 12 * labelSide;
            const w = ctx.measureText(label).width + 6;
            ctx.fillStyle = 'rgba(255,255,255,0.85)';
            ctx.fillRect(lx - w / 2, ly - 8, w, 16);
            ctx.fillStyle = color;
            ctx.fillText(label, lx, ly);
            ctx.textBaseline = 'alphabetic';
        }
    }

    vector(x0, y0, fx, fy, color, name, labelSide = 1) {
        const n = Math.hypot(fx, fy);
        if (n < 0.05) return;
        // fy est vers le haut (physique) -> canvas vers le bas
        this.arrow(x0, y0, fx / n, -fy / n, Animation.arrowLength(n), color,
                   steady(`force-${name}`, () => `${name} ${n.toFixed(1)} N`), labelSide);
    }

    drawForceVectors(cartX, theta, f) {
        const y = this.cartY, L = 2 * this.currentL * this.scale;
        const gx = cartX + (L / 2) * Math.sin(theta), gy = y - (L / 2) * Math.cos(theta);
        const C = FORCE_COLORS;

        // Sur la tige : poids en G, force du pivot au pivot
        this.vector(gx, gy, 0, -f.weightPole, C.weight, 'P', 1);
        this.vector(cartX, y, f.pivot[0], f.pivot[1], C.pivot, 'R', -1);

        // Sur le chariot : poids, reaction du rail (aux roues), moteur, frottement
        this.vector(cartX, y + 4, 0, -f.weightCart, C.weight, 'Pc', 1);

        // N repartie sur les deux roues (N/2 chacune) : le rail pousse les roues vers le haut,
        // fleches dessinees sous le rail et se terminant au contact
        const railY = this.railY;
        const nLen = Animation.arrowLength(f.normal / 2), nDir = f.normal >= 0 ? -1 : 1;
        for (const k of [-1, 1]) {
            this.arrow(cartX + k * this.cartWidth / 3, railY - nDir * nLen, 0, nDir, nLen, C.normal, '');
        }
        const ctx = this.ctx;
        ctx.font = 'bold 12px sans-serif';
        ctx.textAlign = 'left';
        ctx.fillStyle = C.normal;
        ctx.fillText(steady('force-N', () => `N ${Math.abs(f.normal).toFixed(1)} N`), cartX + this.cartWidth / 3 + 8, railY + nLen * 0.6 + 4);

        // Main qui tire le chariot (au-dessus du chariot)
        if (Math.abs(f.hand) > 0.05) {
            this.vector(cartX, y - this.cartHeight / 2 - 4, f.hand, 0, C.hand, 'Main', -1);
        }
        if (Math.abs(f.motor) > 0.05) {
            const s = Math.sign(f.motor);
            this.vector(cartX + s * this.cartWidth / 2, y - 6, f.motor, 0, C.motor, 'F', -s);
        }
        // Frottement au contact roues/rail, oppose au mouvement
        if (Math.abs(f.friction) > 0.05) {
            const s = Math.sign(f.friction);
            this.vector(cartX + s * this.cartWidth / 2, railY - 4, f.friction, 0, C.friction, 'f', s);
        }
    }

    drawForceLegend() {
        const ctx = this.ctx, C = FORCE_COLORS;
        const items = [
            [C.weight, 'P, Pc : poids (tige, chariot)'],
            [C.motor, 'F : force du moteur'],
            [C.hand, 'Main : votre traction sur le chariot'],
            [C.normal, 'N : réaction du rail'],
            [C.pivot, 'R : force du pivot sur la tige'],
            [C.friction, 'f : frottements du rail']
        ];
        ctx.font = '12px sans-serif';
        ctx.textAlign = 'left';
        const x = 10, y0 = this.canvas.height - 12 - (items.length - 1) * 17;
        items.forEach(([c, t], i) => {
            ctx.fillStyle = c;
            ctx.fillRect(x, y0 + i * 17 - 9, 14, 4);
            ctx.fillStyle = '#2c3e50';
            ctx.fillText(t, x + 20, y0 + i * 17 - 3);
        });
        ctx.fillStyle = '#7f8c8d';
        ctx.textAlign = 'right';
        // Barre d'echelle : 10 N
        const bar = 10 * FORCE_SCALE, bx = this.canvas.width - 10 - bar, by = this.canvas.height - 30;
        ctx.strokeStyle = '#2c3e50';
        ctx.lineWidth = 2;
        ctx.beginPath();
        ctx.moveTo(bx, by - 4); ctx.lineTo(bx, by + 4);
        ctx.moveTo(bx, by); ctx.lineTo(bx + bar, by);
        ctx.moveTo(bx + bar, by - 4); ctx.lineTo(bx + bar, by + 4);
        ctx.stroke();
        ctx.fillStyle = '#2c3e50';
        ctx.textAlign = 'center';
        ctx.fillText('10 N', bx + bar / 2, by - 8);
        ctx.fillStyle = '#7f8c8d';
        ctx.textAlign = 'right';
        ctx.fillText('Longueur des flèches proportionnelle à la force', this.canvas.width - 10, this.canvas.height - 10);
        ctx.textAlign = 'left';
    }

    drawRail() {
        const ctx = this.ctx;
        ctx.strokeStyle = '#34495e';
        ctx.lineWidth = 3;
        ctx.beginPath();
        ctx.moveTo(50, this.railY);
        ctx.lineTo(this.canvas.width - 50, this.railY);
        ctx.stroke();

        // Graduations tous les 0.5 m
        ctx.fillStyle = '#95a5a6';
        ctx.font = '11px sans-serif';
        ctx.textAlign = 'center';
        for (let m = -2; m <= 2; m += 0.5) {
            const px = this.centerX + m * this.scale;
            ctx.fillRect(px - 1, this.railY, 2, m % 1 === 0 ? 10 : 6);
            if (m % 1 === 0) ctx.fillText(`${m} m`, px, this.railY + 24);
        }
        ctx.textAlign = 'left';
    }

    drawCart(x) {
        const ctx = this.ctx, w = this.cartWidth, h = this.cartHeight, y = this.cartY;
        const hl = this.dragTarget === 'cart' || (!this.isDragging && this.hitCart(this.mouseX, this.mouseY));

        ctx.fillStyle = hl ? '#2980b9' : '#3498db';
        ctx.fillRect(x - w / 2, y - h / 2, w, h);
        if (hl) {
            ctx.strokeStyle = '#1abc9c';
            ctx.lineWidth = 2;
            ctx.strokeRect(x - w / 2, y - h / 2, w, h);
        }

        ctx.fillStyle = '#2c3e50';
        for (const dx of [-w / 3, w / 3]) {
            ctx.beginPath();
            ctx.arc(x + dx, y + h / 2, this.wheelRadius, 0, 2 * Math.PI);
            ctx.fill();
        }
    }

    // Angle theta : verticale de reference (haut) en pointilles + arc jusqu'a la tige
    drawAngle(cartX, theta) {
        const ctx = this.ctx, y = this.cartY;
        const L = 2 * this.currentL * this.scale;
        const up = -Math.PI / 2;                 // direction "vers le haut" dans le canvas
        const r = Math.max(35, Math.min(90, L * 0.6));

        ctx.save();
        ctx.strokeStyle = '#7f8c8d';
        ctx.lineWidth = 1.5;
        ctx.setLineDash([6, 5]);
        ctx.beginPath();
        ctx.moveTo(cartX, y);
        ctx.lineTo(cartX, y - L - 15);
        ctx.stroke();
        ctx.setLineDash([]);

        if (Math.abs(theta) > 0.003) {
            // Arc entre la verticale et la tige (theta > 0 = sens horaire a l'ecran)
            ctx.strokeStyle = ctx.fillStyle = '#8e44ad';
            ctx.lineWidth = 2;
            ctx.beginPath();
            ctx.arc(cartX, y, r, up, up + theta, theta < 0);
            ctx.stroke();

            // Pointe de fleche au bout de l'arc
            const end = up + theta, s = Math.sign(theta);
            const ex = cartX + r * Math.cos(end), ey = y + r * Math.sin(end);
            const tx = -Math.sin(end) * s, ty = Math.cos(end) * s;   // tangente dans le sens de l'arc
            const nx = Math.cos(end), ny = Math.sin(end);
            ctx.beginPath();
            ctx.moveTo(ex, ey);
            ctx.lineTo(ex - 8 * tx + 4 * nx, ey - 8 * ty + 4 * ny);
            ctx.lineTo(ex - 8 * tx - 4 * nx, ey - 8 * ty - 4 * ny);
            ctx.closePath();
            ctx.fill();
        }

        // Etiquette au milieu de l'arc ; pour un petit angle, de l'autre cote de la verticale
        // pour ne pas chevaucher la tige
        // Hysteresis : le cote ne change que si theta depasse nettement 0, sinon l'etiquette
        // sauterait d'un cote a l'autre a chaque petite oscillation
        if (Math.abs(theta) > 0.05) this.angleLabelSide = -Math.sign(theta);
        const side = this.angleLabelSide || 1;
        const mid = up + (Math.abs(theta) > 0.8 ? theta / 2 : side * 0.45);
        const lr = r + 22 + 30 * Math.abs(Math.cos(mid)); // s'eloigne du chariot quand l'etiquette est laterale
        ctx.fillStyle = '#8e44ad';
        ctx.font = 'bold 13px sans-serif';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'middle';
        ctx.fillText(steady('theta-arc', () => `θ = ${(theta * 180 / Math.PI).toFixed(1)}°`), cartX + lr * Math.cos(mid), y + lr * Math.sin(mid));
        ctx.restore();
    }

    drawPole(cartX) {
        const ctx = this.ctx;
        const bob = this.bobPx();

        const hl = this.dragTarget === 'pendulum' || (!this.isDragging && this.hitBob(this.mouseX, this.mouseY));

        // Tige homogene
        ctx.strokeStyle = hl ? '#c0392b' : '#e74c3c';
        ctx.lineWidth = 9;
        ctx.lineCap = 'round';
        ctx.beginPath();
        ctx.moveTo(cartX, this.cartY);
        ctx.lineTo(bob.x, bob.y);
        ctx.stroke();
        ctx.lineCap = 'butt';

        // Centre de masse (a l du pivot)
        ctx.fillStyle = '#fff';
        ctx.beginPath();
        ctx.arc((cartX + bob.x) / 2, (this.cartY + bob.y) / 2, 3, 0, 2 * Math.PI);
        ctx.fill();

        // Poignee a l'extremite
        ctx.fillStyle = hl ? '#c0392b' : '#e74c3c';
        ctx.beginPath();
        ctx.arc(bob.x, bob.y, 9, 0, 2 * Math.PI);
        ctx.fill();
        if (hl) {
            ctx.strokeStyle = '#1abc9c';
            ctx.lineWidth = 2;
            ctx.stroke();
        }

        ctx.fillStyle = '#2c3e50';
        ctx.beginPath();
        ctx.arc(cartX, this.cartY, 6, 0, 2 * Math.PI);
        ctx.fill();
    }

    drawForce(cartX, F) {
        if (Math.abs(F) < 0.05) return;
        const ctx = this.ctx, s = Math.sign(F);
        const len = Math.min(120, Math.abs(F) * 3);
        const x0 = cartX - s * (this.cartWidth / 2 + len + 4);
        const x1 = cartX - s * (this.cartWidth / 2 + 4);
        const y = this.cartY;
        ctx.strokeStyle = ctx.fillStyle = '#27ae60';
        ctx.lineWidth = 3;
        ctx.beginPath();
        ctx.moveTo(x0, y);
        ctx.lineTo(x1, y);
        ctx.stroke();
        ctx.beginPath();
        ctx.moveTo(x1, y);
        ctx.lineTo(x1 - 9 * s, y - 6);
        ctx.lineTo(x1 - 9 * s, y + 6);
        ctx.closePath();
        ctx.fill();
    }

    drawTarget(x) {
        const ctx = this.ctx;
        ctx.strokeStyle = '#27ae60';
        ctx.lineWidth = 2;
        ctx.setLineDash([5, 5]);
        ctx.beginPath();
        ctx.moveTo(x, this.railY - 120);
        ctx.lineTo(x, this.railY + 10);
        ctx.stroke();
        ctx.setLineDash([]);
    }

    drawValues(x, theta, F) {
        const ctx = this.ctx;
        ctx.fillStyle = '#2c3e50';
        ctx.font = '14px monospace';
        ctx.fillText(steady('hud-x', () => `x = ${x.toFixed(3)} m`), 10, 25);
        ctx.fillText(steady('hud-th', () => `θ = ${theta.toFixed(3)} rad (${(theta * 180 / Math.PI).toFixed(1)}°)`), 10, 45);
        ctx.fillText(steady('hud-F', () => `F = ${F.toFixed(2)} N`), 10, 65);
    }

    drawInstructions() {
        const ctx = this.ctx;
        ctx.fillStyle = '#7f8c8d';
        ctx.font = '12px sans-serif';
        ctx.textAlign = 'right';
        ['Glisser le chariot ou la masse', '← / → : pousser le chariot']
            .forEach((t, i) => ctx.fillText(t, this.canvas.width - 10, 20 + i * 15));
        ctx.textAlign = 'left';
    }
}

// ===== MODULE PLOTS =====
class Plotter {
    constructor(canvasId, title, unit, color, timeWindow = 10) {
        this.canvas = document.getElementById(canvasId);
        this.ctx = this.canvas.getContext('2d');
        this.title = title;
        this.unit = unit;
        this.color = color;
        this.timeWindow = timeWindow;
        this.data = [];
        this.reference = null;
    }

    addPoint(t, y) {
        // Sous-echantillonnage : 1 point toutes les 20 ms suffit pour l'affichage
        const last = this.data[this.data.length - 1];
        if (last && t - last.t < 0.02) return;
        this.data.push({ t, y });
        while (this.data.length && this.data[0].t < t - this.timeWindow) this.data.shift();
    }

    draw() {
        const ctx = this.ctx, W = this.canvas.width, H = this.canvas.height;
        ctx.clearRect(0, 0, W, H);

        const m = { top: 30, right: 20, bottom: 40, left: 60 };
        const pw = W - m.left - m.right, ph = H - m.top - m.bottom;

        const tNow = this.data.length ? this.data[this.data.length - 1].t : 0;
        const tStart = Math.max(0, tNow - this.timeWindow);
        const values = this.data.map(d => d.y);
        if (this.reference !== null) values.push(this.reference);
        // Axe Y stable : bornes arrondies a un pas "rond", elargies tout de suite,
        // resserrees seulement quand la courbe occupe moins de 40 % de l'axe
        const lo = Math.min(0, ...values), hi = Math.max(0, ...values);
        const step = Plotter.niceStep((hi - lo) || 0.01);
        const tMin = Math.floor(lo / step) * step - step / 2, tMax = Math.ceil(hi / step) * step + step / 2;
        if (this.yMin === undefined || lo < this.yMin || hi > this.yMax ||
            (hi - lo) < 0.4 * (this.yMax - this.yMin)) {
            this.yMin = tMin;
            this.yMax = tMax;
        }
        const yMin = this.yMin, yMax = this.yMax;

        const X = (t) => m.left + pw * (t - tStart) / this.timeWindow;
        const Y = (v) => m.top + ph * (yMax - v) / (yMax - yMin);

        // Axes + graduations
        ctx.strokeStyle = '#34495e';
        ctx.lineWidth = 1;
        ctx.beginPath();
        ctx.moveTo(m.left, m.top);
        ctx.lineTo(m.left, m.top + ph);
        ctx.lineTo(m.left + pw, m.top + ph);
        ctx.stroke();

        ctx.fillStyle = '#7f8c8d';
        ctx.font = '12px sans-serif';
        ctx.textAlign = 'right';
        for (let i = 0; i <= 4; i++) {
            const v = yMax - (yMax - yMin) * i / 4;
            ctx.fillText(v.toFixed(2), m.left - 5, Y(v) + 4);
        }
        ctx.textAlign = 'center';
        for (let i = 0; i <= 5; i++) {
            const t = tStart + this.timeWindow * i / 5;
            ctx.fillText(t.toFixed(1), X(t), m.top + ph + 18);
        }
        ctx.fillText('Temps (s)', m.left + pw / 2, H - 5);
        ctx.save();
        ctx.translate(15, m.top + ph / 2);
        ctx.rotate(-Math.PI / 2);
        ctx.fillText(this.unit, 0, 0);
        ctx.restore();

        // Zero et reference
        ctx.strokeStyle = '#ecf0f1';
        ctx.beginPath();
        ctx.moveTo(m.left, Y(0));
        ctx.lineTo(m.left + pw, Y(0));
        ctx.stroke();
        if (this.reference !== null) {
            ctx.strokeStyle = '#27ae60';
            ctx.setLineDash([5, 5]);
            ctx.beginPath();
            ctx.moveTo(m.left, Y(this.reference));
            ctx.lineTo(m.left + pw, Y(this.reference));
            ctx.stroke();
            ctx.setLineDash([]);
        }

        // Courbe
        if (this.data.length > 1) {
            ctx.strokeStyle = this.color;
            ctx.lineWidth = 2;
            ctx.beginPath();
            this.data.forEach((d, i) => (i ? ctx.lineTo(X(d.t), Y(d.y)) : ctx.moveTo(X(d.t), Y(d.y))));
            ctx.stroke();
        }

        ctx.fillStyle = '#2c3e50';
        ctx.font = 'bold 16px sans-serif';
        ctx.fillText(this.title, W / 2, 20);
    }

    static niceStep(range) {
        const raw = range / 4;
        const p = Math.pow(10, Math.floor(Math.log10(raw)));
        const f = raw / p;
        return (f <= 1 ? 1 : f <= 2 ? 2 : f <= 5 ? 5 : 10) * p;
    }

    clear() {
        this.data = [];
        this.yMin = this.yMax = undefined;
    }
}

// ===== MODULE PERFORMANCE =====
class PerformanceMonitor {
    constructor() {
        this.reset();
    }

    reset(mode = 'strict', initialX = 0) {
        this.mode = mode;
        this.initialX = initialX;
        this.stabilizationTime = null;
        this.stableSince = null;
        this.riseTime = null;
        this.maxAngle = 0;
        this.steadyError = null;
        this.fallen = false;
        this.positions = [];
    }

    update(t, x, theta, targetX) {
        const th = Math.abs(theta);
        this.maxAngle = Math.max(this.maxAngle, th);
        if (th > Math.PI / 2) this.fallen = true;

        // Stabilise = |theta| < seuil et x proche de la cible pendant 1 s
        const thetaTol = this.mode === 'strict' ? 0.05 : 0.35;
        const ok = th < thetaTol && Math.abs(x - targetX) < 0.05;
        if (!ok) {
            this.stableSince = null;
            if (this.stabilizationTime !== null && th > 0.35) this.stabilizationTime = null;
        } else if (this.stableSince === null) {
            this.stableSince = t;
        } else if (this.stabilizationTime === null && t - this.stableSince >= 1.0) {
            this.stabilizationTime = this.stableSince;
        }

        if (this.mode === 'project') {
            const span = targetX - this.initialX;
            if (this.riseTime === null && Math.abs(span) > 1e-6 && (x - this.initialX) / span >= 0.9) {
                this.riseTime = t;
            }
            // Erreur statique moyenne sur la derniere seconde
            this.positions.push({ t, x });
            while (this.positions.length && this.positions[0].t < t - 1) this.positions.shift();
            if (t > 2 && Math.abs(span) > 1e-6) {
                const avg = this.positions.reduce((s, p) => s + p.x, 0) / this.positions.length;
                this.steadyError = Math.abs(avg - targetX) / Math.abs(span) * 100;
            }
        }
    }

    getMetrics(force) {
        const status = (v, thr) => (v === null ? null : v <= thr ? 'success' : 'danger');
        const project = this.mode === 'project';

        const m = {
            force,
            stabilityTime: this.stabilizationTime,
            stabilityStatus: status(this.stabilizationTime, 5),
            maxAngle: this.maxAngle,
            angleStatus: status(this.maxAngle, project ? 0.35 : 0.2),
            riseTime: this.riseTime,
            riseStatus: project ? status(this.riseTime, 2) : null,
            steadyError: this.steadyError,
            errorStatus: project ? status(this.steadyError, 2) : null
        };

        if (this.fallen) {
            m.systemStatus = 'Chute détectée';
            m.overallStatus = 'danger';
        } else if (this.stabilizationTime !== null) {
            m.systemStatus = 'Stabilisé';
            const allOk = [m.stabilityStatus, m.angleStatus, m.riseStatus, m.errorStatus]
                .every(s => s === null || s === 'success');
            m.overallStatus = allOk ? 'success' : 'warning';
        } else {
            m.systemStatus = 'En cours';
            m.overallStatus = 'warning';
        }
        return m;
    }
}

// ===== BOUCLE DE COMMANDE NUMERIQUE =====
// Comme sur un vrai banc : le calculateur echantillonne les capteurs toutes les Ts,
// maintient la commande constante entre deux echantillons (bloqueur d'ordre 0),
// et le moteur ne suit la consigne qu'avec un retard du 1er ordre (constante tau).
class ControlLoop {
    constructor(physics, controller, sensors) {
        this.physics = physics;
        this.controller = controller;
        this.sensors = sensors;
        this.controller.plant = physics;
        this.active = true;      // false : controleur coupe, le pendule reste libre
        this.enabled = true;
        this.ts = 0.01;
        this.tau = 0.02;
        this.reset();
    }

    configure({ enabled, noise, tau, ts }) {
        this.enabled = enabled;
        this.sensors.enabled = enabled;
        this.sensors.noiseTheta = noise;
        this.tau = tau;
        this.ts = ts;
    }

    reset() {
        this.controller.reset();
        this.sensors.reset();
        this.clock = 0;
        this.command = 0;
        this.motorForce = 0;
    }

    step(dt, xRef, push = 0, integrator = 'rk4', hand = null) {
        const p = this.physics;

        // Chariot tenu a la main : le moteur est coupe, la physique continue
        if (hand !== null) {
            this.command = 0;
            this.motorForce = 0;
            p.step(push, integrator, hand);
            return;
        }

        if (!this.active) {
            this.command = 0;
            this.motorForce = 0;
            p.step(push, integrator);
            return;
        }

        if (!this.enabled) {
            this.command = this.controller.compute(p.state, xRef, dt);
            p.step(this.command + push, integrator);
            return;
        }

        this.clock -= dt;
        if (this.clock <= 0) {
            const Ts = Math.max(this.ts, dt);
            this.clock += Ts;
            const measured = this.sensors.read(p.state, Ts);
            this.command = clamp(this.controller.compute(measured, xRef, Ts), -p.fMax, p.fMax);
        }

        this.motorForce += (this.command - this.motorForce) * (this.tau > 0 ? 1 - Math.exp(-dt / this.tau) : 1);
        p.step(this.motorForce + push, integrator);
    }
}

// ===== APPLICATION PRINCIPALE =====
class InvertedPendulumApp {
    constructor() {
        this.physics = new Physics();
        this.controller = new PIDController();
        this.ui = new UI();
        this.animation = new Animation('animation-canvas');
        this.anglePlot = new Plotter('angle-plot', 'Angle θ(t)', 'θ (rad)', '#e74c3c');
        this.positionPlot = new Plotter('position-plot', 'Position x(t)', 'x (m)', '#3498db');
        this.performance = new PerformanceMonitor();
        this.bench = new ControlLoop(this.physics, this.controller, new Sensors());

        this.isRunning = false;
        this.animationId = null;
        this.lastTime = 0;
        this.accumulator = 0;
        this.pendingImpulseAt = null;
        this.keys = { left: false, right: false };
        this.keyForce = 15; // N
        this.handTargetX = null; // position visee par la main quand on tire le chariot

        this.bindUI();
        this.bindInteraction();
        this.bindKeyboard();

        // Gains de depart calcules pour les parametres affiches
        this.autoTune();
        this.reset();
        this.play();
    }

    bindUI() {
        this.ui.onPlay = () => this.start();
        this.ui.onPause = () => this.pause();
        this.ui.onReset = () => this.reset();
        this.ui.onImpulse = () => this.applyImpulse();
        this.ui.onAutoTune = () => this.autoTune();
        this.ui.onPresetChange = () => this.applyPreset();
        this.ui.onScenarioChange = () => this.loadScenario();
        this.ui.onModeChange = () => this.reset();
    }

    bindInteraction() {
        // Pendant un glissement la simulation est figee, puis repart a la relache
        // Glisser le chariot : la souris le tire comme une main (ressort amorti vers la souris).
        // La simulation continue, donc le pendule reagit aux accelerations du chariot.
        this.animation.onCartDrag = (dx) => {
            const lim = this.physics.railLimit;
            if (this.handTargetX === null) this.handTargetX = this.physics.x;
            this.handTargetX = clamp(this.handTargetX + dx, -lim, lim);
        };
        this.animation.onPendulumDrag = (theta) => {
            this.physics.setState(this.physics.x, 0, theta, 0);
        };
        this.animation.onDragEnd = () => {
            this.handTargetX = null;
            this.bench.reset();
            this.measuring = false;
        };
    }

    bindKeyboard() {
        const set = (e, v) => {
            if (e.target && ['INPUT', 'SELECT'].includes(e.target.tagName)) return;
            if (e.key === 'ArrowLeft') { this.keys.left = v; e.preventDefault(); }
            if (e.key === 'ArrowRight') { this.keys.right = v; e.preventDefault(); }
        };
        window.addEventListener('keydown', (e) => set(e, true));
        window.addEventListener('keyup', (e) => set(e, false));
    }

    get targetX() {
        return this.ui.mode === 'project' ? 0.2 : 0;
    }

    // « Lancer » : active le controleur (redressement puis stabilisation)
    start() {
        this.bench.active = true;
        this.play();
    }

    play() {
        if (this.isRunning) return;
        this.isRunning = true;
        this.lastTime = performance.now();
        this.accumulator = 0;
        this.animationId = requestAnimationFrame((t) => this.loop(t));
    }

    pause() {
        this.isRunning = false;
        if (this.animationId) cancelAnimationFrame(this.animationId);
        this.animationId = null;
    }

    reset() {
        const scenario = this.ui.scenarioSelect.value;
        // Redressement : le pendule pend au repos et le controleur attend « Lancer ».
        // Impulsion : on part a l'equilibre haut et on perturbe apres 1 s.
        const initial = {
            swingup: [0, 0, Math.PI, 0],
            impulse: [0, 0, 0, 0]
        }[scenario] || [0, 0, 0.1, 0];

        this.physics.updateParams(this.ui.getSystemParams());
        this.physics.reset(initial);
        this.bench.reset();
        this.bench.active = scenario !== 'swingup';
        this.measuring = false;
        this.forceDisplay = null;
        this.performance.reset(this.ui.mode, 0);
        this.performance.t0 = 0;
        this.pendingImpulseAt = scenario === 'impulse' ? 1.0 : null;

        this.anglePlot.clear();
        this.positionPlot.clear();
        this.positionPlot.reference = this.ui.mode === 'project' ? this.targetX : null;

        this.render();
    }

    applyImpulse() {
        this.physics.applyImpulse(1.0); // 1 N.s
    }

    autoTune() {
        const p = this.ui.getSystemParams();
        const speed = PIDController.presetSpeed(this.ui.pidPreset.value);
        this.ui.setGains(PIDController.autoTune(p, speed));
        if (this.ui.pidPreset.value === 'custom') this.ui.pidPreset.value = 'nominal';
    }

    applyPreset() {
        if (this.ui.pidPreset.value !== 'custom') this.autoTune();
    }

    loadScenario() {
        if (this.ui.scenarioSelect.value === 'position') this.ui.setMode('project');
        else this.ui.setMode('strict');
        this.reset();
        this.play();
    }

    loop(now) {
        if (!this.isRunning) return;

        // Pas fixe, synchronise sur le temps reel (max 0.1 s rattrape par image)
        const sim = this.ui.getSimulationParams();
        this.physics.dt = sim.dt;
        this.accumulator += Math.min((now - this.lastTime) / 1000, 0.1);
        this.lastTime = now;

        // Seule la masse tenue a la main fige la simulation ; le chariot tire reste simule
        if (this.animation.dragTarget !== 'pendulum') {
            const r = this.ui.getRealismParams();
            this.physics.updateParams({ ...this.ui.getSystemParams(), coulomb: r.enabled ? r.fs : 0 });
            this.controller.setGains(this.ui.getGains());
            this.bench.configure(r);
            while (this.accumulator >= sim.dt) {
                this.update(sim);
                this.accumulator -= sim.dt;
            }
        } else {
            this.accumulator = 0;
        }

        this.render();
        this.animationId = requestAnimationFrame((t) => this.loop(t));
    }

    update(sim) {
        const p = this.physics;

        if (this.pendingImpulseAt !== null && p.t >= this.pendingImpulseAt) {
            this.applyImpulse();
            this.pendingImpulseAt = null;
            // Les mesures demarrent a la perturbation
            this.performance.reset(this.ui.mode, p.x);
            this.performance.t0 = p.t;
        }

        let push = 0;
        if (this.keys.left) push -= this.keyForce;
        if (this.keys.right) push += this.keyForce;

        this.bench.step(sim.dt, this.targetX, push, sim.integrator, this.handForce());

        // Les indicateurs mesurent la stabilisation : ils demarrent quand le pendule est rattrape en haut
        const balancing = this.bench.active && this.controller.phase === 'balance';
        if (balancing && !this.measuring) {
            this.measuring = true;
            this.performance.reset(this.ui.mode, p.x);
            this.performance.t0 = p.t;
        } else if (!balancing) {
            this.measuring = false;
        }
        if (this.measuring) {
            this.performance.update(p.t - this.performance.t0, p.x, p.theta, this.targetX);
        }
        this.anglePlot.addPoint(p.t, p.theta);
        this.positionPlot.addPoint(p.t, p.x);
        this.smoothForces(p.forces(p.lastForce), sim.dt);
    }

    // Main : ressort-amortisseur entre le chariot et la souris (frequence ~4 Hz, bien amorti)
    handForce() {
        if (this.handTargetX === null) return null;
        const p = this.physics, M = p.mc + p.mp;
        const w = 2 * Math.PI * 4;
        return M * w * w * (this.handTargetX - p.x) - 2 * 0.8 * M * w * p.x_dot;
    }

    // Filtre passe-bas (constante 0.1 s) sur les forces affichees : le bruit des capteurs fait
    // varier la commande a chaque echantillon, ce qui ferait clignoter fleches et valeurs.
    smoothForces(f, dt) {
        const d = this.forceDisplay;
        if (!d) {
            this.forceDisplay = { ...f, pivot: [...f.pivot] };
            return;
        }
        const a = 1 - Math.exp(-dt / 0.1);
        for (const k of ['motor', 'hand', 'weightPole', 'weightCart', 'normal', 'friction']) d[k] += a * (f[k] - d[k]);
        d.pivot[0] += a * (f.pivot[0] - d.pivot[0]);
        d.pivot[1] += a * (f.pivot[1] - d.pivot[1]);
    }

    render() {
        const p = this.physics;
        const showForces = document.getElementById('show-forces').checked;
        if (!this.forceDisplay || this.animation.dragTarget === 'pendulum') {
            this.forceDisplay = null;
            this.smoothForces(p.forces(p.lastForce), 0);
        }
        const motor = this.forceDisplay.motor;
        this.animation.draw(p.x, p.theta, p.l, motor, this.ui.mode === 'project' ? this.targetX : null,
                            showForces ? this.forceDisplay : null);
        this.anglePlot.draw();
        this.positionPlot.draw();
        const m = this.performance.getMetrics(motor);
        if (!this.measuring) m.angleStatus = null;
        if (!this.bench.active) {
            m.systemStatus = 'Au repos';
            m.overallStatus = null;
        } else if (this.controller.phase === 'swing') {
            m.systemStatus = 'Redressement';
            m.overallStatus = 'warning';
        }
        this.ui.updateIndicators(m);
    }
}

// ===== INITIALISATION =====
if (typeof document !== 'undefined') {
    document.addEventListener('DOMContentLoaded', () => {
        window.pendulumApp = new InvertedPendulumApp();
    });
}

if (typeof module !== 'undefined') {
    module.exports = { Physics, PIDController, Sensors, ControlLoop, wrapAngle };
}
