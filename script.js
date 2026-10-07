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
// Tous les textes changent dans la meme image (meme tranche de 250 ms) : ils restent coherents entre eux.
const steady = (key, make) => {
    const slot = Math.floor(Date.now() / TEXT_REFRESH_MS);
    const c = textCache.get(key);
    if (c && c.slot === slot) return c.s;
    const s = make();
    textCache.set(key, { slot, s });
    return s;
};

const wrapAngle = (a) => {
    const twoPi = 2 * Math.PI;
    a = (a + Math.PI) % twoPi;
    if (a < 0) a += twoPi;
    return a - Math.PI;
};

// Geometrie du pendule (longueur totale de la tige L = 2l) :
//  - 'bob' : masse ponctuelle au bout d'une tige legere -> centre de masse a L, J = mp L^2
//  - 'rod' : tige homogene                               -> centre de masse a L/2, J = mp L^2 / 3
function poleGeometry(p) {
    const L = 2 * p.l;
    return p.poleType === 'rod'
        ? { L, lc: L / 2, J: p.mp * L * L / 3 }
        : { L, lc: L, J: p.mp * L * L };
}

// ===== MODULE PHYSICS =====
// Banc de type laboratoire : chariot entraine par un moteur a courant continu (pignon/cremaillere),
// rail fini avec butees amorties, pendule = tige homogene. Tous les parametres sont en SI.
class Physics {
    constructor(params = {}) {
        this.mc = 1.0;        // masse chariot + inertie ramenee du moteur (kg)
        this.mp = 0.25;       // masse de la tige (kg)
        this.l = 0.3;         // demi-longueur de la tige (m) : tige de 2l = 0.6 m
        this.poleType = 'bob';  // 'bob' : masse au bout ; 'rod' : tige homogene
        this.b = 0.5;         // frottement visqueux du chariot (N.s/m)
        this.c = 0.002;       // frottement du pivot (N.m.s)
        this.g = 9.81;        // gravite (m/s2)
        this.coulomb = 0;     // frottement sec du chariot (N)
        this.railLimit = 0.6; // course du centre du chariot : +-0.6 m (contact avec les butees)

        // Moteur DC + reducteur + pignon : F = alpha V - beta x' (force contre-electromotrice)
        this.motorModel = true;
        this.alpha = 1.72;    // N/V
        this.beta = 7.7;      // N.s/m
        this.vMax = 10;       // tension max de l'amplificateur (V)

        // Butees : ressort + amortisseur raides (caoutchouc)
        this.kStop = 2e4;     // N/m
        this.cStop = 120;     // N.s/m

        Object.assign(this, params);

        this.state = [0, 0, 0.1, 0]; // [x, x_dot, theta, theta_dot]
        this.t = 0;
        this.dt = 0.005;
        this.input = { cmd: 0, cart: 0, tip: [0, 0] };
        this.lastForce = 0;
    }

    updateParams(params) {
        Object.assign(this, params);
    }

    // Force max a l'arret (pour l'affichage et le mode ideal)
    get fMax() { return this.alpha * this.vMax; }

    get geom() { return poleGeometry(this); }

    // Force du moteur selon la commande et la vitesse du chariot
    motorForce(x_dot) {
        if (this.motorModel) {
            const V = clamp(this.input.cmd, -this.vMax, this.vMax);
            return this.alpha * V - this.beta * x_dot;
        }
        return clamp(this.input.cmd, -this.fMax, this.fMax); // source de force ideale
    }

    // Force des butees (nulle tant que le chariot ne les touche pas)
    stopForce(x, x_dot) {
        const pen = Math.abs(x) - this.railLimit;
        if (pen <= 0) return 0;
        const s = Math.sign(x);
        // L'amortisseur ne tire jamais le chariot vers la butee
        return Math.min(0, -this.kStop * pen - this.cStop * x_dot * s) * s;
    }

    // Equations de Lagrange avec forces generalisees :
    //  (mc+mp) x'' + mp l cos th th'' = Qx + mp l th'^2 sin th
    //  cos th x'' + Le th''           = g sin th - c/(mp l) th' + Qth/(mp l)
    // Qx = moteur + frottements + butees + main (chariot) + main (bout de tige)
    // Qth = couple de la main au bout de la tige (point a 2l du pivot)
    derivatives(state) {
        const [x, x_dot, theta, theta_dot] = state;
        const S = Math.sin(theta), C = Math.cos(theta);
        const [fx, fy] = this.input.tip;

        const Qx = this.motorForce(x_dot) - this.b * x_dot - this.coulomb * Math.tanh(x_dot / 0.01)
                 + this.stopForce(x, x_dot) + this.input.cart + fx;
        const { L, lc, J } = this.geom;
        const ml = this.mp * lc;
        const Qth = L * (fx * C - fy * S);

        // Equation en theta divisee par mp lc : longueur equivalente Le = J / (mp lc)
        const a11 = this.mc + this.mp, a12 = ml * C;
        const a21 = C, a22 = J / ml;
        const v1 = Qx + ml * theta_dot * theta_dot * S;
        const v2 = this.g * S - (this.c / ml) * theta_dot + Qth / ml;

        const det = a11 * a22 - a12 * a21; // = (mc+mp) Le - mp l cos^2 > 0
        return [x_dot, (a22 * v1 - a12 * v2) / det, theta_dot, (-a21 * v1 + a11 * v2) / det];
    }

    rk4Step(h) {
        const s = this.state;
        const k1 = this.derivatives(s);
        const k2 = this.derivatives(s.map((v, i) => v + h * k1[i] / 2));
        const k3 = this.derivatives(s.map((v, i) => v + h * k2[i] / 2));
        const k4 = this.derivatives(s.map((v, i) => v + h * k3[i]));
        for (let i = 0; i < 4; i++) s[i] += h * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6;
    }

    eulerStep(h) {
        const d = this.derivatives(this.state);
        for (let i = 0; i < 4; i++) this.state[i] += h * d[i];
    }

    // cmd : tension (V) avec le modele moteur, sinon force (N)
    // ext.cart : force exterieure sur le chariot (N) ; ext.tip : [fx, fy] au bout de la tige (N)
    step(cmd, integrator = 'rk4', ext = {}) {
        this.input = { cmd, cart: ext.cart || 0, tip: ext.tip || [0, 0] };

        // Pres des butees, sous-pas pour integrer le choc raide sans instabilite
        const nearStop = Math.abs(this.state[0]) > this.railLimit - 0.05;
        const n = nearStop ? 8 : 1, h = this.dt / n;
        for (let i = 0; i < n; i++) {
            if (integrator === 'rk4') this.rk4Step(h);
            else this.eulerStep(h);
        }

        this.state[2] = wrapAngle(this.state[2]);
        this.lastForce = this.motorForce(this.state[1]);
        this.t += this.dt;
    }

    reset(initialState = [0, 0, 0.1, 0]) {
        this.state = [...initialState];
        this.t = 0;
        this.input = { cmd: 0, cart: 0, tip: [0, 0] };
        this.lastForce = 0;
    }

    // Impulsion horizontale (N.s) sur le chariot : partagee avec le pendule via les equations
    applyImpulse(impulse) {
        this.state[1] += impulse / this.mc;
    }

    // Bilan des forces a l'instant courant (N), repere x vers la droite, y vers le haut.
    // Force du pivot sur la tige par Newton au centre de masse G (main au bout comprise) :
    //   R + mp g + F_tip = mp aG
    forces() {
        const [x, x_dot, theta, theta_dot] = this.state;
        const [, x_ddot, , theta_ddot] = this.derivatives(this.state);
        const S = Math.sin(theta), C = Math.cos(theta);
        const [fx, fy] = this.input.tip;
        const lc = this.geom.lc;
        const aGx = x_ddot + lc * (C * theta_ddot - S * theta_dot ** 2);
        const aGy = -lc * (S * theta_ddot + C * theta_dot ** 2);
        const Rx = this.mp * aGx - fx;
        const Ry = this.mp * aGy + this.mp * this.g - fy;
        return {
            motor: this.motorForce(x_dot),
            hand: this.input.cart,
            tip: [fx, fy],
            stop: this.stopForce(x, x_dot),
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
        this.swingGain = 6;      // gain du pompage d'energie
        this.swingAccel = 0.5;   // acceleration max du chariot en redressement (en g)
        this.swingKx = 1;        // rappel du chariot vers le centre pendant le redressement
        this.swingKv = 3;
        this.swingMargin = 0.2;  // vise un peu plus que l'energie du sommet pour l'atteindre
        this.swingWall = 0.5;    // intensite du mur de fin de course (en g)
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
        const { mc, mp, g } = this.plant;
        const { lc, J } = poleGeometry(this.plant);     // J : inertie autour du pivot
        const E = 0.5 * J * theta_dot ** 2 + mp * g * lc * (Math.cos(theta) - 1);
        const E0 = 2 * mp * g * lc;

        // Au repos parfait en bas, sign(0) = 0 : on donne un premier coup
        const dir = Math.abs(theta_dot) > 1e-3 ? Math.sign(theta_dot * Math.cos(theta)) : 1;
        const aMax = this.swingAccel * g;
        let a = clamp(this.swingGain * g * (E / E0 - this.swingMargin) * dir, -aMax, aMax);

        // Rappel doux du chariot vers le centre pour rester sur le rail
        a += -this.swingKx * (x - xRef) - this.swingKv * x_dot;

        // Mur progressif : au-dela de la moitie de la course, ramene fermement le chariot
        // vers le centre pour ne pas taper les butees
        const R = this.plant.railLimit || 1;
        const over = clamp((Math.abs(x) - 0.5 * R) / (0.5 * R), 0, 1);
        a -= Math.sign(x) * this.swingWall * g * over * over;

        return (mc + mp) * a;
    }

    reset() {
        this.integral = 0;
        this.phase = 'swing';
    }

    // Auto-reglage par placement de poles (Ackermann) sur le modele linearise.
    // speed > 1 = reponse plus rapide (et forces plus grandes).
    static autoTune(p, speed = 1) {
        const { mc, mp, b, c, g } = p;
        const { lc: l, J } = poleGeometry(p);
        const M = mc + mp, Le = J / (mp * l);
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
const SLIDERS = ['kp', 'ki', 'kd', 'kx', 'kv', 'mc', 'mp', 'l', 'b', 'c', 'g', 'vmax', 'rail', 'dt',
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
        this.$('pole-type').addEventListener('change', () => this.onPoleChange());
    }

    refreshValue(name) {
        const v = parseFloat(this.$(`${name}-slider`).value);
        const unit = { dt: ' ms', ts: ' ms', tau: ' ms', noise: ' mrad', fs: ' N', vmax: ' V', rail: ' m' }[name];
        this.$(`${name}-value`).textContent = unit ? `${+v.toFixed(2)}${unit}` : v.toFixed(name === 'c' ? 3 : 2);
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
    onPoleChange() {}
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
        return { mc: v('mc'), mp: v('mp'), l: v('l') / 2, b: v('b'), c: v('c'), g: v('g'),
                 poleType: this.$('pole-type').value,
                 vMax: v('vmax'), railLimit: v('rail') };
    }

    getSimulationParams() {
        return {
            dt: parseFloat(this.$('dt-slider').value) / 1000, // ms -> s
            integrator: this.integratorSelect.value
        };
    }
}

// ===== MODULE ANIMATION AVEC INTERACTIVITE =====
// Dimensions reelles du banc (m)
const CART_W = 0.20, CART_H = 0.07, WHEEL_R = 0.018;

const FORCE_SCALE = 12;    // px par newton
const FORCE_MAX_PX = 200;  // au-dela, la fleche est tronquee (la valeur reste affichee)
const FORCE_COLORS = {
    weight: '#d35400',   // poids
    motor: '#27ae60',    // moteur
    hand: '#2c3e50',     // main (glisser le chariot ou la tige)
    stop: '#c0392b',     // butee du rail
    tension: '#8e44ad',  // force de la tige sur la masse
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

        this.onDragStart = null;
        this.onDragEnd = null;
        this.dragWorld = null;

        this.currentX = 0;
        this.currentTheta = 0;
        this.currentL = 0.3;
        this.poleType = 'bob';  // fixe par l'application
        this.railLimit = 0.6;   // course du centre du chariot (m), fixee par l'application

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

    // Pixels du canvas -> metres (origine : milieu du rail, hauteur du pivot)
    toWorld(px, py) {
        return { x: (px - this.centerX) / this.scale, y: (this.cartY - py) / this.scale };
    }

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
        this.dragWorld = this.toWorld(x, y);
        if (this.onDragStart) this.onDragStart(this.dragTarget, this.dragWorld);
        this.canvas.style.cursor = 'grabbing';
        this.mouseX = x; this.mouseY = y;
    }

    move({ x, y }) {
        if (!this.isDragging) {
            this.canvas.style.cursor = (this.hitBob(x, y) || this.hitCart(x, y)) ? 'grab' : 'default';
            this.mouseX = x; this.mouseY = y;
            return;
        }
        this.dragWorld = this.toWorld(x, y);
        this.mouseX = x; this.mouseY = y;
    }

    end() {
        if (this.isDragging && this.onDragEnd) this.onDragEnd();
        this.isDragging = false;
        this.dragTarget = null;
        this.dragWorld = null;
        this.canvas.style.cursor = 'default';
    }

    draw(x, theta, l, force, targetX = null, forces = null, grabOffset = 0) {
        this.currentX = x;
        this.currentTheta = theta;
        this.currentL = l;
        this.layout(l);

        const ctx = this.ctx;
        ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);

        this.drawRail();
        if (targetX !== null) this.drawTarget(this.centerX + targetX * this.scale);

        const cartX = this.cartPx();
        if (!forces) this.drawForce(cartX, force);
        this.drawCart(cartX);
        this.drawAngle(cartX, theta);
        this.drawPole(cartX);
        this.drawHand(grabOffset);
        if (forces) {
            this.drawForceVectors(cartX, theta, forces);
            this.drawForceLegend();
        }
        this.drawValues(x, theta, force);
        this.drawInstructions();
    }

    // Echelle commune a tout le dessin : le rail complet tient en largeur, la tige en hauteur
    layout(l) {
        const W = this.canvas.width, H = this.canvas.height;
        const half = this.railLimit + CART_W / 2;           // position des butees
        this.scale = Math.min((W - 170) / (2 * half), (this.cartY - 75) / (2 * l));
        this.cartWidth = CART_W * this.scale;
        this.cartHeight = Math.max(14, CART_H * this.scale);
        this.wheelRadius = Math.max(5, WHEEL_R * this.scale);
        this.railY = this.cartY + this.cartHeight / 2 + this.wheelRadius;
        this.stopPx = half * this.scale;
    }

    // Longueur proportionnelle a la force : FORCE_SCALE px par newton
    static arrowLength(F) {
        return Math.min(FORCE_MAX_PX, FORCE_SCALE * Math.abs(F));
    }

    // Fleche partant de (x0, y0) dans la direction (ux, uy) (repere canvas), de longueur len
    arrow(x0, y0, ux, uy, len, color, label, labelSide = 1, dashed = false) {
        if (len < 1.5) return;
        const ctx = this.ctx;
        const x1 = x0 + ux * len, y1 = y0 + uy * len;
        ctx.strokeStyle = ctx.fillStyle = color;
        ctx.lineWidth = dashed ? 1.5 : 2.5;
        if (dashed) ctx.setLineDash([5, 4]);
        ctx.beginPath();
        ctx.moveTo(x0, y0);
        ctx.lineTo(x1, y1);
        ctx.stroke();
        ctx.setLineDash([]);
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

    vector(x0, y0, fx, fy, color, name, labelSide = 1, dashed = false) {
        const n = Math.hypot(fx, fy);
        if (n < 0.05) return;
        // fy est vers le haut (physique) -> canvas vers le bas
        this.arrow(x0, y0, fx / n, -fy / n, Animation.arrowLength(n), color,
                   steady(`force-${name}`, () => `${name} ${n.toFixed(1)} N`), labelSide, dashed);
    }

    drawForceVectors(cartX, theta, f) {
        const y = this.cartY, L = 2 * this.currentL * this.scale;
        const bob = this.poleType !== 'rod';
        const S = Math.sin(theta), Cs = Math.cos(theta);
        // Point d'application du poids : la masse au bout, ou le milieu de la tige homogene
        const gl = bob ? L : L / 2;
        const gx = cartX + gl * S, gy = y - gl * Cs;
        const C = FORCE_COLORS;

        // Poids du pendule et sa decomposition (le long de la tige / perpendiculaire)
        const P = f.weightPole;
        const Ppar = -P * Cs;                         // composante sur u = (sin, cos), axe pivot -> masse
        // (inutile quand la tige est presque verticale : P∥ se confond avec P).
        // Pour la masse au bout, la decomposition est dans l'encart zoome.
        if (!bob && Math.abs(S) > 0.08) {
            this.vector(gx, gy, Ppar * S, Ppar * Cs, C.weight, 'P∥', -1, true);
            this.vector(gx, gy, -Ppar * S, -P - Ppar * Cs, C.weight, 'P⊥', 1, true);
        }
        this.vector(gx, gy, 0, -P, C.weight, 'P', 1);

        if (bob) {
            // Masse au bout : la tige (legere) tire ou pousse la masse le long de son axe : T
            this.vector(gx, gy, f.pivot[0], f.pivot[1], C.tension, 'T', -1);
        }
        // Force du chariot sur la tige, au pivot
        this.vector(cartX, y, f.pivot[0], f.pivot[1], C.pivot, 'R', -1);

        if (bob) this.drawBobInset(theta, f);

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

        // Main au bout de la tige
        const tip = this.bobPx();
        this.vector(tip.x, tip.y, f.tip[0], f.tip[1], C.hand, 'Main', 1);
        // Butee : pousse le chariot vers l'interieur au contact
        if (Math.abs(f.stop) > 0.05) {
            const s = Math.sign(f.stop);
            this.vector(cartX - s * this.cartWidth / 2, y, f.stop, 0, C.stop, 'Butée', 1);
        }
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

    // Encart zoome : bilan des forces sur la masse (diagramme du corps isole), a sa propre echelle
    drawBobInset(theta, f) {
        const ctx = this.ctx, C = FORCE_COLORS;
        const size = 210, x0 = this.canvas.width - size - 10, y0 = 50;
        const cx = x0 + size / 2, cy = y0 + size / 2 + 8;
        const S = Math.sin(theta), Cs = Math.cos(theta);
        const P = f.weightPole, Ppar = -P * Cs;
        const vecs = [
            { v: [0, -P], c: C.weight, n: 'P', side: 1 },
            { v: [f.pivot[0], f.pivot[1]], c: C.tension, n: 'T', side: -1 },
            { v: [Ppar * S, Ppar * Cs], c: C.weight, n: 'P∥', side: -1, dashed: true },
            { v: [-Ppar * S, -P - Ppar * Cs], c: C.weight, n: 'P⊥', side: 1, dashed: true },
            { v: f.tip, c: C.hand, n: 'Main', side: 1 }
        ];
        // Echelle de l'encart : la plus grande force occupe ~75 px
        const maxF = Math.max(0.5, ...vecs.map(o => Math.hypot(o.v[0], o.v[1])));
        const k = Math.min(40, 75 / maxF);

        ctx.save();
        ctx.fillStyle = 'rgba(255,255,255,0.92)';
        ctx.strokeStyle = '#bdc3c7';
        ctx.lineWidth = 1;
        ctx.fillRect(x0, y0, size, size);
        ctx.strokeRect(x0, y0, size, size);
        ctx.fillStyle = '#2c3e50';
        ctx.font = 'bold 12px sans-serif';
        ctx.textAlign = 'left';
        ctx.fillText('Forces sur la masse (zoom)', x0 + 8, y0 + 16);

        // Direction de la tige (vers le pivot) et masse
        ctx.strokeStyle = '#7f8c8d';
        ctx.lineWidth = 3;
        ctx.beginPath();
        ctx.moveTo(cx, cy);
        ctx.lineTo(cx - 70 * S, cy + 70 * Cs);
        ctx.stroke();
        ctx.fillStyle = '#e74c3c';
        ctx.beginPath(); ctx.arc(cx, cy, 10, 0, 2 * Math.PI); ctx.fill();
        ctx.restore();

        for (const o of vecs) {
            const n = Math.hypot(o.v[0], o.v[1]);
            if (n < 0.02 || (o.dashed && Math.abs(S) < 0.08)) continue;
            this.arrow(cx, cy, o.v[0] / n, -o.v[1] / n, n * k, o.c,
                       steady(`inset-${o.n}`, () => `${o.n} ${n.toFixed(2)} N`), o.side, o.dashed);
        }
    }

    // Barre d'echelle des forces (la legende des couleurs est dans la page, sous la scene)
    drawForceLegend() {
        const ctx = this.ctx;
        ctx.font = '12px sans-serif';
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
        const ctx = this.ctx, cx = this.centerX, y = this.railY, sp = this.stopPx;

        // Courroie (moteur a gauche, poulie de renvoi a droite), accrochee au chariot
        const pulleyR = 9, motorX = cx - sp - 34, idlerX = cx + sp + 22, beltY = y + 10;
        ctx.strokeStyle = '#7f8c8d';
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        ctx.moveTo(motorX, beltY - pulleyR); ctx.lineTo(idlerX, beltY - pulleyR);
        ctx.moveTo(motorX, beltY + pulleyR); ctx.lineTo(idlerX, beltY + pulleyR);
        ctx.stroke();
        ctx.fillStyle = '#95a5a6';
        ctx.beginPath(); ctx.arc(idlerX, beltY, pulleyR, 0, 2 * Math.PI); ctx.fill();

        // Moteur a courant continu
        ctx.fillStyle = '#2c3e50';
        ctx.fillRect(motorX - 30, beltY - 16, 30, 32);
        ctx.fillStyle = '#95a5a6';
        ctx.beginPath(); ctx.arc(motorX, beltY, pulleyR, 0, 2 * Math.PI); ctx.fill();
        ctx.fillStyle = '#fff';
        ctx.font = 'bold 12px sans-serif';
        ctx.textAlign = 'center';
        ctx.fillText('M', motorX - 15, beltY + 4);

        // Rail
        ctx.fillStyle = '#34495e';
        ctx.fillRect(cx - sp - 10, y, 2 * sp + 20, 4);

        // Butees (bloc + tampon caoutchouc)
        for (const k of [-1, 1]) {
            const ex = cx + k * sp;
            ctx.fillStyle = '#2c3e50';
            ctx.fillRect(k < 0 ? ex - 12 : ex + 4, y - 26, 8, 30);
            ctx.fillStyle = '#e67e22';
            ctx.fillRect(k < 0 ? ex - 4 : ex, y - 20, 4, 18);
        }

        // Graduations tous les 10 cm, etiquettes tous les 20 cm
        ctx.fillStyle = '#95a5a6';
        ctx.font = '11px sans-serif';
        const n = Math.floor(this.railLimit * 10 + 1e-6);
        for (let i = -n; i <= n; i++) {
            const px = cx + i * 0.1 * this.scale;
            ctx.fillRect(px - 0.5, y + 4, 1, i % 2 === 0 ? 8 : 5);
            if (i % 2 === 0) ctx.fillText(`${(i / 10).toFixed(1)}`, px, y + 34);
        }
        ctx.fillText('x (m)', cx + sp + 30, y + 34);
        ctx.textAlign = 'left';
    }

    drawCart(x) {
        const ctx = this.ctx, w = this.cartWidth, h = this.cartHeight, y = this.cartY;
        const hl = this.dragTarget === 'cart' || (!this.isDragging && this.hitCart(this.mouseX, this.mouseY));

        ctx.fillStyle = hl ? '#2980b9' : '#3498db';
        ctx.fillRect(x - w / 2, y - h / 2, w, h);
        ctx.strokeStyle = hl ? '#1abc9c' : '#21618c';
        ctx.lineWidth = hl ? 2 : 1;
        ctx.strokeRect(x - w / 2, y - h / 2, w, h);

        // Roues qui tournent : angle = distance parcourue / rayon
        const r = this.wheelRadius, rot = this.currentX / WHEEL_R;
        for (const dx of [-w / 3, w / 3]) {
            const wx = x + dx, wy = y + h / 2;
            ctx.fillStyle = '#2c3e50';
            ctx.beginPath(); ctx.arc(wx, wy, r, 0, 2 * Math.PI); ctx.fill();
            ctx.strokeStyle = '#bdc3c7';
            ctx.lineWidth = 1.5;
            ctx.beginPath();
            for (let k = 0; k < 3; k++) {
                const a = rot + k * Math.PI / 3;
                ctx.moveTo(wx - r * 0.8 * Math.cos(a), wy - r * 0.8 * Math.sin(a));
                ctx.lineTo(wx + r * 0.8 * Math.cos(a), wy + r * 0.8 * Math.sin(a));
            }
            ctx.stroke();
        }
    }

    // Ligne elastique entre le point saisi et la souris (la "main")
    drawHand(grabOffset) {
        const w = this.dragWorld;
        if (!w) return;
        const ctx = this.ctx;
        const from = this.dragTarget === 'cart'
            ? { x: this.cartPx() + grabOffset * this.scale, y: this.cartY }
            : this.bobPx();
        const to = { x: this.centerX + w.x * this.scale, y: this.cartY - w.y * this.scale };
        ctx.save();
        ctx.strokeStyle = '#2c3e50';
        ctx.lineWidth = 1.5;
        ctx.setLineDash([4, 4]);
        ctx.beginPath(); ctx.moveTo(from.x, from.y); ctx.lineTo(to.x, to.y); ctx.stroke();
        ctx.setLineDash([]);
        ctx.fillStyle = 'rgba(44,62,80,0.25)';
        ctx.beginPath(); ctx.arc(to.x, to.y, 9, 0, 2 * Math.PI); ctx.fill();
        ctx.restore();
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
        // Etiquette gardee dans le canvas
        const lx = clamp(cartX + lr * Math.cos(mid), 50, this.canvas.width - 50);
        ctx.fillText(steady('theta-arc', () => `θ = ${(theta * 180 / Math.PI).toFixed(1)}°`), lx, y + lr * Math.sin(mid));
        ctx.restore();
    }

    drawPole(cartX) {
        const ctx = this.ctx;
        const bob = this.bobPx();

        const hl = this.dragTarget === 'pendulum' || (!this.isDragging && this.hitBob(this.mouseX, this.mouseY));

        const massive = this.poleType === 'rod';
        // Tige : epaisse et rouge si elle porte la masse, fine et grise si elle est legere
        ctx.strokeStyle = massive ? (hl ? '#c0392b' : '#e74c3c') : '#7f8c8d';
        ctx.lineWidth = massive ? Math.max(4, 0.02 * this.scale) : 3;
        ctx.lineCap = 'round';
        ctx.beginPath();
        ctx.moveTo(cartX, this.cartY);
        ctx.lineTo(bob.x, bob.y);
        ctx.stroke();
        ctx.lineCap = 'butt';

        if (massive) {
            // Centre de masse de la tige homogene (a mi-longueur)
            ctx.fillStyle = '#fff';
            ctx.beginPath();
            ctx.arc((cartX + bob.x) / 2, (this.cartY + bob.y) / 2, 3, 0, 2 * Math.PI);
            ctx.fill();
        }

        // Extremite : poignee (tige homogene) ou masse ponctuelle (boule de 3 cm de rayon)
        const r = massive ? 9 : Math.max(10, 0.03 * this.scale);
        ctx.fillStyle = hl ? '#c0392b' : '#e74c3c';
        ctx.beginPath();
        ctx.arc(bob.x, bob.y, r, 0, 2 * Math.PI);
        ctx.fill();
        if (!massive) {
            ctx.fillStyle = '#fff';
            ctx.beginPath();
            ctx.arc(bob.x, bob.y, 2.5, 0, 2 * Math.PI);
            ctx.fill();
        }
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
        if (this.voltage !== undefined) {
            ctx.fillText(steady('hud-V', () => `V = ${this.voltage.toFixed(1)} V`), 10, 85);
        }
    }

    drawInstructions() {
        const ctx = this.ctx;
        ctx.fillStyle = '#7f8c8d';
        ctx.font = '12px sans-serif';
        ctx.textAlign = 'right';
        ['Saisir le chariot ou le bout de la tige', '← / → : pousser le chariot']
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
        this.ts = 0.005;   // periode d'echantillonnage (s)
        this.tau = 0.005;  // constante de temps amplificateur + moteur (s)
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
        this.voltage = 0;
    }

    // ext : forces exterieures (clavier / main) ; held : l'utilisateur tient le systeme
    step(dt, xRef, ext = {}, integrator = 'rk4', held = false) {
        const p = this.physics;

        // Amplificateur coupe (au repos ou tenu a la main) : le moteur n'agit plus
        if (held || !this.active) {
            this.command = this.voltage = 0;
            p.motorModel = false;
            p.step(0, integrator, ext);
            return;
        }

        // Mode ideal : mesure parfaite, commande continue, source de force parfaite
        if (!this.enabled) {
            p.motorModel = false;
            this.command = this.controller.compute(p.state, xRef, dt);
            p.step(this.command, integrator, ext);
            return;
        }

        // Banc reel : calculateur echantillonne a Ts. Il calcule une force voulue F*,
        // puis la tension V = (F* + beta x'_mesure) / alpha qui compense la force
        // contre-electromotrice du moteur, saturee a +-Vmax.
        p.motorModel = true;
        this.clock -= dt;
        if (this.clock <= 0) {
            const Ts = Math.max(this.ts, dt);
            this.clock += Ts;
            const measured = this.sensors.read(p.state, Ts);
            const Fwanted = this.controller.compute(measured, xRef, Ts);
            this.command = clamp((Fwanted + p.beta * measured[1]) / p.alpha, -p.vMax, p.vMax);
        }

        // Amplificateur + inductance du moteur : retard du 1er ordre sur la tension
        this.voltage += (this.command - this.voltage) * (this.tau > 0 ? 1 - Math.exp(-dt / this.tau) : 1);
        p.step(this.voltage, integrator, ext);
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
        this.keyForce = 8; // N, poussee de la main au clavier
        this.grabOffset = 0;

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
        // Autre modele de pendule : nouveaux gains et depart de zero
        this.ui.onPoleChange = () => { this.autoTune(); this.reset(); };
        this.ui.onModeChange = () => this.reset();
    }

    // La souris agit comme une main : un ressort amorti relie le point saisi a la souris.
    // La simulation ne s'arrete jamais ; le moteur est coupe tant qu'on tient le systeme.
    bindInteraction() {
        this.animation.onDragStart = (target, w) => {
            // Decalage entre le point saisi et le centre du chariot
            this.grabOffset = target === 'cart' ? w.x - this.physics.x : 0;
        };
        this.animation.onDragEnd = () => {
            this.bench.reset();
            this.measuring = false;
        };
    }

    // Forces exterieures de l'utilisateur (clavier + main)
    externalForces() {
        const ext = { cart: 0, tip: [0, 0] };
        if (this.keys.left) ext.cart -= this.keyForce;
        if (this.keys.right) ext.cart += this.keyForce;

        const a = this.animation, w = a.dragWorld, p = this.physics;
        if (!w) return ext;
        if (a.dragTarget === 'cart') {
            // Main sur le chariot : ressort horizontal ~3 Hz, force limitee a 60 N
            const M = p.mc + p.mp, om = 2 * Math.PI * 3;
            const F = M * om * om * (w.x - this.grabOffset - p.x) - 2 * 0.8 * M * om * p.x_dot;
            ext.cart += clamp(F, -60, 60);
        } else if (a.dragTarget === 'pendulum') {
            // Main au bout de la tige : ressort 2D vers la souris, force limitee a 30 N
            const L = 2 * p.l, S = Math.sin(p.theta), C = Math.cos(p.theta);
            const tip = [p.x + L * S, L * C];
            const v = [p.x_dot + L * C * p.theta_dot, -L * S * p.theta_dot];
            let fx = 150 * (w.x - tip[0]) - 6 * v[0];
            let fy = 150 * (w.y - tip[1]) - 6 * v[1];
            const n = Math.hypot(fx, fy);
            if (n > 30) { fx *= 30 / n; fy *= 30 / n; }
            ext.tip = [fx, fy];
        }
        return ext;
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

        const r = this.ui.getRealismParams();
        this.physics.updateParams({ ...this.ui.getSystemParams(), coulomb: r.enabled ? r.fs : 0 });
        this.controller.setGains(this.ui.getGains());
        this.bench.configure(r);
        while (this.accumulator >= sim.dt) {
            this.update(sim);
            this.accumulator -= sim.dt;
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

        this.bench.step(sim.dt, this.targetX, this.externalForces(), sim.integrator, this.animation.isDragging);

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
        this.smoothForces(p.forces(), sim.dt);
    }

    // Filtre passe-bas (constante 0.1 s) sur les forces affichees : le bruit des capteurs fait
    // varier la commande a chaque echantillon, ce qui ferait clignoter fleches et valeurs.
    smoothForces(f, dt) {
        const d = this.forceDisplay;
        if (!d) {
            this.forceDisplay = { ...f, pivot: [...f.pivot], tip: [...f.tip] };
            return;
        }
        const a = 1 - Math.exp(-dt / 0.1);
        for (const k of ['motor', 'hand', 'stop', 'weightPole', 'weightCart', 'normal', 'friction']) d[k] += a * (f[k] - d[k]);
        for (const k of ['pivot', 'tip']) for (const i of [0, 1]) d[k][i] += a * (f[k][i] - d[k][i]);
    }

    render() {
        const p = this.physics;
        const showForces = document.getElementById('show-forces').checked;
        document.getElementById('force-legend').hidden = !showForces;
        if (!this.forceDisplay) this.smoothForces(p.forces(), 0);
        const motor = this.forceDisplay.motor;
        this.animation.railLimit = p.railLimit;
        this.animation.poleType = p.poleType;
        this.animation.voltage = p.motorModel ? this.bench.voltage : undefined;
        this.animation.draw(p.x, p.theta, p.l, motor, this.ui.mode === 'project' ? this.targetX : null,
                            showForces ? this.forceDisplay : null, this.grabOffset);
        this.anglePlot.draw();
        this.positionPlot.draw();
        const m = this.performance.getMetrics(motor);
        if (!this.measuring) m.angleStatus = null;
        if (this.animation.isDragging) {
            m.systemStatus = 'Tenu à la main';
            m.overallStatus = null;
        } else if (!this.bench.active) {
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
