// =====================================================================
// Pendule inverse sur chariot - simulation temps reel
// Convention : x vers la droite, theta = 0 pendule vertical VERS LE HAUT,
// theta > 0 quand le pendule penche vers la droite.
// =====================================================================

const clamp = (v, lo, hi) => Math.max(lo, Math.min(hi, v));

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
        Object.assign(this, params);

        this.state = [0, 0, 0.1, 0]; // [x, x_dot, theta, theta_dot]
        this.t = 0;
        this.dt = 0.005;
        this.lastForce = 0;
    }

    updateParams(params) {
        Object.assign(this, params);
    }

    // Equations de Lagrange (masse ponctuelle a la distance l) :
    //  (mc+mp) x'' + mp l cos(th) th'' = F - b x' + mp l th'^2 sin(th)
    //  cos(th) x'' + l th''           = g sin(th) - c/(mp l) th'
    derivatives(state, F) {
        const [, x_dot, theta, theta_dot] = state;
        const S = Math.sin(theta);
        const C = Math.cos(theta);

        const a11 = this.mc + this.mp;
        const a12 = this.mp * this.l * C;
        const a21 = C;
        const a22 = this.l;

        const v1 = F - this.b * x_dot + this.mp * this.l * theta_dot * theta_dot * S;
        const v2 = this.g * S - (this.c / (this.mp * this.l)) * theta_dot;

        // det = l (mc + mp sin^2) > 0 toujours
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

    step(force, integrator = 'rk4') {
        const F = clamp(force, -this.fMax, this.fMax);
        this.lastForce = F;

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
        this.fallAngle = Math.PI / 2; // au-dela le controleur abandonne
    }

    setGains(gains) {
        Object.assign(this, gains);
    }

    compute(state, xRef, dt) {
        const [x, x_dot, theta, theta_dot] = state;
        const th = wrapAngle(theta);

        if (Math.abs(th) > this.fallAngle) {
            this.integral = 0;
            return 0;
        }

        this.integral = clamp(this.integral + th * dt, -this.maxIntegral, this.maxIntegral);

        return this.kp * th + this.ki * this.integral + this.kd * theta_dot
             + this.kx * (x - xRef) + this.kv * x_dot;
    }

    reset() {
        this.integral = 0;
    }

    // Auto-reglage par placement de poles (Ackermann) sur le modele linearise.
    // speed > 1 = reponse plus rapide (et forces plus grandes).
    static autoTune(p, speed = 1) {
        const { mc, mp, l, b, c, g } = p;
        const cp = c / (mp * l);
        const A = [
            [0, 1, 0, 0],
            [0, -b / mc, -mp * g / mc, mp * cp / mc],
            [0, 0, 0, 1],
            [0, b / (l * mc), (mc + mp) * g / (l * mc), -(mc + mp) * cp / (l * mc)]
        ];
        const B = [0, 1 / mc, 0, -1 / (l * mc)];
        const mul = (X, Y) => X.map(r => Y[0].map((_, j) => r.reduce((s, v, k) => s + v * Y[k][j], 0)));
        const mv = (X, v) => X.map(r => r.reduce((s, a, k) => s + a * v[k], 0));
        const I = [0, 1, 2, 3].map(i => [0, 1, 2, 3].map(j => (i === j ? 1 : 0)));

        // 2 poles lents (position) + 2 poles rapides (angle) lies a la pulsation propre
        const w0 = Math.sqrt(Math.abs(g) * (mc + mp) / (mc * l)) || 1;
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
        const M = [0, 1, 2, 3].map(i => [...cols[i], i === 3 ? 1 : 0]);
        for (let i = 0; i < 4; i++) {
            let piv = i;
            for (let r = i + 1; r < 4; r++) if (Math.abs(M[r][i]) > Math.abs(M[piv][i])) piv = r;
            [M[i], M[piv]] = [M[piv], M[i]];
            for (let r = 0; r < 4; r++) {
                if (r === i) continue;
                const f = M[r][i] / M[i][i];
                for (let k = i; k < 5; k++) M[r][k] -= f * M[i][k];
            }
        }
        const y = M.map((r, i) => r[4] / r[i]);
        const K = [0, 1, 2, 3].map(j => y.reduce((s, v, k) => s + v * phi[k][j], 0));

        // F = -K s
        return { kx: -K[0], kv: -K[1], kp: -K[2], ki: 0, kd: -K[3] };
    }

    static presetSpeed(type) {
        return { soft: 0.7, nominal: 1.0, aggressive: 1.4 }[type] || 1.0;
    }
}

// ===== MODULE UI =====
const SLIDERS = ['kp', 'ki', 'kd', 'kx', 'kv', 'mc', 'mp', 'l', 'b', 'c', 'g', 'fmax', 'dt'];

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
        this.$(`${name}-value`).textContent = name === 'dt' ? `${v} ms` : v.toFixed(2);
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
        el.className = 'indicator' + (status ? ` ${status}` : '');
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
class Animation {
    constructor(canvasId) {
        this.canvas = document.getElementById(canvasId);
        this.ctx = this.canvas.getContext('2d');

        this.scale = 150; // pixels par metre
        this.centerX = this.canvas.width / 2;
        this.cartY = this.canvas.height * 0.68; // le pendule monte au-dessus du chariot

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

    bobPx() {
        const L = this.currentL * this.scale;
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

    draw(x, theta, l, force, targetX = null) {
        this.currentX = x;
        this.currentTheta = theta;
        this.currentL = l;

        const ctx = this.ctx;
        ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);

        this.drawRail();
        if (targetX !== null) this.drawTarget(this.centerX + targetX * this.scale);

        const cartX = this.cartPx();
        this.drawForce(cartX, force);
        this.drawCart(cartX);
        this.drawPole(cartX);
        this.drawValues(x, theta, force);
        this.drawInstructions();
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

    drawPole(cartX) {
        const ctx = this.ctx;
        const bob = this.bobPx();

        ctx.strokeStyle = '#e74c3c';
        ctx.lineWidth = 4;
        ctx.beginPath();
        ctx.moveTo(cartX, this.cartY);
        ctx.lineTo(bob.x, bob.y);
        ctx.stroke();

        const hl = this.dragTarget === 'pendulum' || (!this.isDragging && this.hitBob(this.mouseX, this.mouseY));
        ctx.fillStyle = hl ? '#c0392b' : '#e74c3c';
        ctx.beginPath();
        ctx.arc(bob.x, bob.y, 12, 0, 2 * Math.PI);
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
        ctx.fillText(`x = ${x.toFixed(3)} m`, 10, 25);
        ctx.fillText(`θ = ${theta.toFixed(3)} rad (${(theta * 180 / Math.PI).toFixed(1)}°)`, 10, 45);
        ctx.fillText(`F = ${F.toFixed(2)} N`, 10, 65);
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
        let yMin = Math.min(0, ...values), yMax = Math.max(0, ...values);
        const pad = Math.max((yMax - yMin) * 0.1, 0.01);
        yMin -= pad; yMax += pad;

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

    clear() {
        this.data = [];
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

        this.isRunning = false;
        this.animationId = null;
        this.lastTime = 0;
        this.accumulator = 0;
        this.pendingImpulseAt = null;
        this.keys = { left: false, right: false };
        this.keyForce = 15; // N

        this.bindUI();
        this.bindInteraction();
        this.bindKeyboard();

        // Gains de depart calcules pour les parametres affiches
        this.autoTune();
        this.reset();
        this.play();
    }

    bindUI() {
        this.ui.onPlay = () => this.play();
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
        this.animation.onCartDrag = (dx) => {
            const lim = this.physics.railLimit;
            this.physics.setState(clamp(this.physics.x + dx, -lim, lim), 0, this.physics.theta, 0);
        };
        this.animation.onPendulumDrag = (theta) => {
            this.physics.setState(this.physics.x, 0, clamp(theta, -Math.PI / 2, Math.PI / 2), 0);
        };
        this.animation.onDragEnd = () => {
            this.controller.reset();
            this.performance.reset(this.ui.mode, this.physics.x);
            this.performance.t0 = this.physics.t;
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
        // Impulsion : on part a l'equilibre et on perturbe apres 1 s
        const initial = scenario === 'impulse' ? [0, 0, 0, 0] : [0, 0, 0.1, 0];

        this.physics.updateParams(this.ui.getSystemParams());
        this.physics.reset(initial);
        this.controller.reset();
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

        if (!this.animation.isDragging) {
            this.physics.updateParams(this.ui.getSystemParams());
            this.controller.setGains(this.ui.getGains());
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

        let F = this.controller.compute(p.state, this.targetX, sim.dt);
        if (this.keys.left) F -= this.keyForce;
        if (this.keys.right) F += this.keyForce;

        p.step(F, sim.integrator);

        this.performance.update(p.t - (this.performance.t0 || 0), p.x, p.theta, this.targetX);
        this.anglePlot.addPoint(p.t, p.theta);
        this.positionPlot.addPoint(p.t, p.x);
    }

    render() {
        const p = this.physics;
        this.animation.draw(p.x, p.theta, p.l, p.lastForce, this.ui.mode === 'project' ? this.targetX : null);
        this.anglePlot.draw();
        this.positionPlot.draw();
        this.ui.updateIndicators(this.performance.getMetrics(p.lastForce));
    }
}

// ===== INITIALISATION =====
if (typeof document !== 'undefined') {
    document.addEventListener('DOMContentLoaded', () => {
        window.pendulumApp = new InvertedPendulumApp();
    });
}

if (typeof module !== 'undefined') {
    module.exports = { Physics, PIDController, wrapAngle };
}
