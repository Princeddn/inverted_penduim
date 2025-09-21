// ===== MODULE PHYSICS =====
class Physics {
    constructor(params = {}) {
        // Parametres physiques du systeme
        this.mc = params.mc || 1.0;     // masse chariot (kg)
        this.mp = params.mp || 0.3;     // masse pendule (kg)
        this.l = params.l || 0.5;       // longueur pendule (m)
        this.b = params.b || 0.05;      // coefficient frottement chariot (N·s/m)
        this.c = params.c || 0.01;      // coefficient frottement pendule (N·m·s)
        this.g = params.g || 9.81;      // gravite (m/s2)
        this.fMax = params.fMax || 50;  // force max (N)

        // Etat du systeme [x, x_dot, theta, theta_dot]
        this.state = [0, 0, 0.1, 0]; // Commence pres de la verticale

        // Temps
        this.t = 0;
        this.dt = 0.01; // pas d'integration (s)
    }

    // Met a jour les parametres physiques
    updateParams(params) {
        Object.assign(this, params);
    }

    // Normalise l'angle entre -PI et PI (utilitaire du modele de reference)
    wrapAngle(angle) {
        const pi = Math.PI;
        const twoPi = 2 * pi;
        angle = (angle + pi) % twoPi;
        if (angle < 0) angle += twoPi;
        return angle - pi;
    }

    // Calcule les derivees du systeme - Equations du modele de reference
    derivatives(state, force) {
        const [x, x_dot, theta, theta_dot] = state;
        const F = Math.max(-this.fMax, Math.min(this.fMax, force)); // saturation

        const S = Math.sin(theta);
        const C = Math.cos(theta);
        const totalM = this.mc + this.mp;

        // Frottement pendule (ajoute pour plus de realisme)
        const c = this.c || 0.01; // coefficient frottement pendule

        // Matrice implicite 2x2 sur [x¨, θ¨]^T (methode de reference)
        // (mc+mp) * x¨ + b x˙ + mp l (θ¨ cosθ - θ˙² sinθ) = F
        // l θ¨ + (c/(mp l)) θ˙ + g sinθ = - x¨ cosθ

        const a11 = totalM;
        const a12 = this.mp * this.l * C;
        const a21 = C;
        const a22 = this.l;

        const v1 = F - this.b * x_dot - this.mp * this.l * theta_dot * theta_dot * S;
        const v2 = -(c / (this.mp * this.l)) * theta_dot - this.g * S;

        // Resolution 2x2 : det = (mc+mp)*l - mp*l*C*C = l(mc + mp*(1 - C^2))
        const det = a11 * a22 - a12 * a21;
        const inv11 = a22 / det;
        const inv12 = -a12 / det;
        const inv21 = -a21 / det;
        const inv22 = a11 / det;

        const x_ddot = inv11 * v1 + inv12 * v2;
        const theta_ddot = inv21 * v1 + inv22 * v2;

        return [x_dot, x_ddot, theta_dot, theta_ddot];
    }

    // Integration Runge-Kutta 4eme ordre
    rk4Step(force) {
        const k1 = this.derivatives(this.state, force);

        const state_k2 = this.state.map((s, i) => s + this.dt * k1[i] / 2);
        const k2 = this.derivatives(state_k2, force);

        const state_k3 = this.state.map((s, i) => s + this.dt * k2[i] / 2);
        const k3 = this.derivatives(state_k3, force);

        const state_k4 = this.state.map((s, i) => s + this.dt * k3[i]);
        const k4 = this.derivatives(state_k4, force);

        // Mise a jour de l'etat
        for (let i = 0; i < 4; i++) {
            this.state[i] += this.dt * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]) / 6;
        }

        // Normalise l'angle theta (indice 2) pour eviter les problemes numeriques
        this.state[2] = this.wrapAngle(this.state[2]);

        this.t += this.dt;
    }

    // Integration Euler
    eulerStep(force) {
        const derivatives = this.derivatives(this.state, force);
        for (let i = 0; i < 4; i++) {
            this.state[i] += this.dt * derivatives[i];
        }
        this.t += this.dt;
    }

    // Avance la simulation d'un pas
    step(force, integrator = 'rk4') {
        if (integrator === 'rk4') {
            this.rk4Step(force);
        } else {
            this.eulerStep(force);
        }

        // Limite la position du chariot avec amortissement intelligent
        const maxPosition = 1.1; // metres (reduit pour plus de securite)
        const softBoundary = 0.9; // zone d'amortissement

        if (this.state[0] > maxPosition) {
            this.state[0] = maxPosition;
            this.state[1] = -Math.abs(this.state[1]) * 0.3; // Rebond tres amorti
        } else if (this.state[0] < -maxPosition) {
            this.state[0] = -maxPosition;
            this.state[1] = Math.abs(this.state[1]) * 0.3; // Rebond tres amorti
        }

        // Zone d'amortissement progressive plus aggressive
        if (this.state[0] > softBoundary) {
            const factor = (this.state[0] - softBoundary) / (maxPosition - softBoundary);
            this.state[1] *= (1 - factor * 0.6); // Amortissement plus fort
        } else if (this.state[0] < -softBoundary) {
            const factor = (-this.state[0] - softBoundary) / (maxPosition - softBoundary);
            this.state[1] *= (1 - factor * 0.6); // Amortissement plus fort
        }

        // Limitation de la vitesse angulaire pour eviter le spinning
        const maxAngularVelocity = 10; // rad/s
        if (Math.abs(this.state[3]) > maxAngularVelocity) {
            this.state[3] = Math.sign(this.state[3]) * maxAngularVelocity;
        }
    }

    // Remet a zero la simulation
    reset(initialState = [0, 0, 0.1, 0]) {
        this.state = [...initialState];
        this.t = 0;
    }

    // Applique une impulsion
    applyImpulse(impulse) {
        // Impulsion = changement de momentum
        // Pour le chariot: Δp = mc * Δv = impulse
        this.state[1] += impulse / this.mc;
    }

    // Verifie si le pendule est tombe (avec detection de boucle infinie)
    hasFallen() {
        // Detection de problemes graves
        const isOutOfBounds = Math.abs(this.state[0]) > 1.8;
        const isSpinning = Math.abs(this.state[2]) > 6*Math.PI; // Plus de 3 tours
        const hasExtremeVelocity = Math.abs(this.state[1]) > 10 || Math.abs(this.state[3]) > 20;

        return isOutOfBounds || isSpinning || hasExtremeVelocity;
    }

    // Permet de modifier directement l'etat (pour l'interactivite)
    setState(x, x_dot, theta, theta_dot) {
        this.state = [x, x_dot, theta, theta_dot];
    }

    // Getters pour l'etat actuel
    get x() { return this.state[0]; }
    get x_dot() { return this.state[1]; }
    get theta() { return this.state[2]; }
    get theta_dot() { return this.state[3]; }
}

// ===== MODULE CONTROLLER =====
class PIDController {
    constructor(kp = 80, ki = 15, kd = 12) {
        this.kp = kp;
        this.ki = ki;
        this.kd = kd;

        this.integral = 0;
        this.lastError = 0;
        this.lastTime = 0;
    }

    // Met a jour les gains PID
    setGains(kp, ki, kd) {
        this.kp = kp;
        this.ki = ki;
        this.kd = kd;
    }

    // Calcule la commande PID
    compute(setpoint, measurement, dt) {
        const error = setpoint - measurement;

        // Terme integral
        this.integral += error * dt;

        // Terme derive
        const derivative = (error - this.lastError) / dt;

        // Commande PID
        const output = this.kp * error + this.ki * this.integral + this.kd * derivative;

        this.lastError = error;

        return output;
    }

    // Remet a zero l'integrateur
    reset() {
        this.integral = 0;
        this.lastError = 0;
    }

    // Auto-reglage Ziegler-Nichols (approximation)
    autoTune() {
        // Approximation simple basee sur les caracteristiques du systeme
        const ku = 25; // Gain critique approximatif
        const tu = 0.5; // Periode critique approximative

        this.kp = 0.6 * ku;
        this.ki = 2 * this.kp / tu;
        this.kd = this.kp * tu / 8;

        return { kp: this.kp, ki: this.ki, kd: this.kd };
    }

    // Presets de reglage
    static getPreset(type) {
        const presets = {
            underdamped: { kp: 45, ki: 10, kd: 5 },
            critical: { kp: 50, ki: 12, kd: 6 },
            overdamped: { kp: 35, ki: 8, kd: 4 }
        };
        return presets[type] || presets.critical;
    }
}

// ===== MODULE UI =====
class UI {
    constructor() {
        this.isPlaying = false;
        this.mode = 'strict'; // 'strict' ou 'project'
        this.initializeElements();
        this.setupEventListeners();
    }

    initializeElements() {
        // Controles principaux
        this.playBtn = document.getElementById('play-btn');
        this.pauseBtn = document.getElementById('pause-btn');
        this.resetBtn = document.getElementById('reset-btn');
        this.impulseBtn = document.getElementById('impulse-btn');

        // Mode
        this.strictModeBtn = document.getElementById('strict-mode');
        this.projectModeBtn = document.getElementById('project-mode');

        // Sliders PID
        this.kpSlider = document.getElementById('kp-slider');
        this.kiSlider = document.getElementById('ki-slider');
        this.kdSlider = document.getElementById('kd-slider');
        this.pidPreset = document.getElementById('pid-preset');
        this.autoTuneBtn = document.getElementById('auto-tune');

        // Sliders systeme
        this.mcSlider = document.getElementById('mc-slider');
        this.mpSlider = document.getElementById('mp-slider');
        this.lSlider = document.getElementById('l-slider');
        this.bSlider = document.getElementById('b-slider');
        this.gSlider = document.getElementById('g-slider');
        this.fmaxSlider = document.getElementById('fmax-slider');

        // Simulation
        this.dtSlider = document.getElementById('dt-slider');
        this.integratorSelect = document.getElementById('integrator');
        this.scenarioSelect = document.getElementById('scenario');

        // Indicateurs
        this.stabilityIndicator = document.getElementById('stability-indicator');
        this.angleIndicator = document.getElementById('angle-indicator');
        this.riseIndicator = document.getElementById('rise-indicator');
        this.errorIndicator = document.getElementById('error-indicator');
        this.statusIndicator = document.getElementById('status-indicator');

        // Valeurs affichees
        this.setupValueDisplays();
    }

    setupValueDisplays() {
        const sliders = [
            'kp', 'ki', 'kd', 'mc', 'mp', 'l', 'b', 'g', 'fmax', 'dt'
        ];

        sliders.forEach(name => {
            const slider = document.getElementById(`${name}-slider`);
            const display = document.getElementById(`${name}-value`);

            slider.addEventListener('input', () => {
                let value = parseFloat(slider.value);
                if (name === 'dt') {
                    display.textContent = value.toString();
                } else {
                    display.textContent = value.toFixed(2);
                }
            });
        });
    }

    setupEventListeners() {
        // Controles principaux
        this.playBtn.addEventListener('click', () => this.onPlay());
        this.pauseBtn.addEventListener('click', () => this.onPause());
        this.resetBtn.addEventListener('click', () => this.onReset());
        this.impulseBtn.addEventListener('click', () => this.onImpulse());

        // Mode
        this.strictModeBtn.addEventListener('click', () => this.setMode('strict'));
        this.projectModeBtn.addEventListener('click', () => this.setMode('project'));

        // Auto-tune
        this.autoTuneBtn.addEventListener('click', () => this.onAutoTune());

        // Presets PID
        this.pidPreset.addEventListener('change', () => this.onPresetChange());

        // Scenarios
        this.scenarioSelect.addEventListener('change', () => this.onScenarioChange());
    }

    setMode(mode) {
        this.mode = mode;
        this.strictModeBtn.classList.toggle('active', mode === 'strict');
        this.projectModeBtn.classList.toggle('active', mode === 'project');
    }

    // Callbacks (a implementer par l'application principale)
    onPlay() { console.log('Play'); }
    onPause() { console.log('Pause'); }
    onReset() { console.log('Reset'); }
    onImpulse() { console.log('Impulse'); }
    onAutoTune() { console.log('Auto-tune'); }
    onPresetChange() { console.log('Preset change'); }
    onScenarioChange() { console.log('Scenario change'); }

    // Met a jour les indicateurs de performance
    updateIndicators(metrics) {
        this.updateIndicator('stability', metrics.stabilityTime, metrics.stabilityStatus);
        this.updateIndicator('angle', metrics.maxAngle, metrics.angleStatus);
        this.updateIndicator('rise', metrics.riseTime, metrics.riseStatus);
        this.updateIndicator('error', metrics.steadyError, metrics.errorStatus);
        this.updateIndicator('status', metrics.systemStatus, metrics.overallStatus);
    }

    updateIndicator(type, value, status) {
        const indicator = document.getElementById(`${type}-indicator`);
        const valueElement = document.getElementById(`${type === 'status' ? 'system-status' : type + '-time'}`);

        if (type === 'status') {
            valueElement.textContent = value;
        } else if (type === 'angle') {
            document.getElementById('max-angle').textContent = `${value.toFixed(3)} rad`;
        } else if (type === 'error') {
            document.getElementById('steady-error').textContent = `${value.toFixed(1)}%`;
        } else {
            valueElement.textContent = value !== null ? `${value.toFixed(2)} s` : '-- s';
        }

        // Met a jour les classes CSS
        indicator.className = 'indicator';
        if (status === 'success') indicator.classList.add('success');
        else if (status === 'warning') indicator.classList.add('warning');
        else if (status === 'danger') indicator.classList.add('danger');
    }

    // Getters pour les parametres
    getPIDParams() {
        return {
            kp: parseFloat(this.kpSlider.value),
            ki: parseFloat(this.kiSlider.value),
            kd: parseFloat(this.kdSlider.value)
        };
    }

    getSystemParams() {
        return {
            mc: parseFloat(this.mcSlider.value),
            mp: parseFloat(this.mpSlider.value),
            l: parseFloat(this.lSlider.value),
            b: parseFloat(this.bSlider.value),
            g: parseFloat(this.gSlider.value),
            fMax: parseFloat(this.fmaxSlider.value)
        };
    }

    getSimulationParams() {
        return {
            dt: parseFloat(this.dtSlider.value) / 1000, // ms vers s
            integrator: this.integratorSelect.value
        };
    }

    // Setters pour mettre a jour l'interface
    setPIDParams(params) {
        this.kpSlider.value = params.kp;
        this.kiSlider.value = params.ki;
        this.kdSlider.value = params.kd;

        // Declenche les evenements pour mettre a jour l'affichage
        this.kpSlider.dispatchEvent(new Event('input'));
        this.kiSlider.dispatchEvent(new Event('input'));
        this.kdSlider.dispatchEvent(new Event('input'));

        this.pidPreset.value = 'custom';
    }
}

// ===== MODULE ANIMATION AVEC INTERACTIVITE =====
class Animation {
    constructor(canvasId) {
        this.canvas = document.getElementById(canvasId);
        this.ctx = this.canvas.getContext('2d');

        // Parametres d'affichage
        this.scale = 150; // pixels par metre (reduit pour garder le pendule visible)
        this.centerX = this.canvas.width / 2;
        this.centerY = this.canvas.height / 2; // Centre verticalement

        // Parametres visuels
        this.cartWidth = 60;
        this.cartHeight = 30;
        this.poleWidth = 4;
        this.railY = this.centerY + this.cartHeight / 2;

        // Variables d'interactivite
        this.isDragging = false;
        this.dragTarget = null; // 'cart' ou 'pendulum'
        this.lastMouseX = 0;
        this.lastMouseY = 0;
        this.mouseForce = 0;

        // Callbacks pour l'interactivite
        this.onCartDrag = null;
        this.onPendulumDrag = null;
        this.onMouseForce = null;

        this.setupInteraction();
    }

    setupInteraction() {
        // Gestion des evenements souris
        this.canvas.addEventListener('mousedown', (e) => this.onMouseDown(e));
        this.canvas.addEventListener('mousemove', (e) => this.onMouseMove(e));
        this.canvas.addEventListener('mouseup', (e) => this.onMouseUp(e));
        this.canvas.addEventListener('mouseleave', (e) => this.onMouseUp(e));

        // Gestion tactile pour mobile
        this.canvas.addEventListener('touchstart', (e) => this.onTouchStart(e));
        this.canvas.addEventListener('touchmove', (e) => this.onTouchMove(e));
        this.canvas.addEventListener('touchend', (e) => this.onTouchEnd(e));

        // Style du curseur
        this.canvas.style.cursor = 'grab';
    }

    getMousePos(e) {
        const rect = this.canvas.getBoundingClientRect();
        return {
            x: e.clientX - rect.left,
            y: e.clientY - rect.top
        };
    }

    getTouchPos(e) {
        const rect = this.canvas.getBoundingClientRect();
        return {
            x: e.touches[0].clientX - rect.left,
            y: e.touches[0].clientY - rect.top
        };
    }

    isPointInCart(x, y, cartX, cartY) {
        return x >= cartX - this.cartWidth / 2 &&
               x <= cartX + this.cartWidth / 2 &&
               y >= cartY - this.cartHeight / 2 &&
               y <= cartY + this.cartHeight / 2;
    }

    isPointInPendulum(x, y, pendulumX, pendulumY) {
        const radius = 15; // Rayon de detection autour de la masse
        const dx = x - pendulumX;
        const dy = y - pendulumY;
        return dx * dx + dy * dy <= radius * radius;
    }

    onMouseDown(e) {
        const pos = this.getMousePos(e);
        this.startInteraction(pos.x, pos.y);
    }

    onTouchStart(e) {
        e.preventDefault();
        const pos = this.getTouchPos(e);
        this.startInteraction(pos.x, pos.y);
    }

    startInteraction(x, y) {
        // Calcule les positions actuelles
        const cartX = this.centerX + this.currentX * this.scale;
        const cartY = this.centerY;

        const poleLength = this.currentL * this.scale;
        const poleEndX = cartX + poleLength * Math.sin(this.currentTheta);
        const poleEndY = cartY - poleLength * Math.cos(this.currentTheta);

        // Detecte quel element est clique
        if (this.isPointInPendulum(x, y, poleEndX, poleEndY)) {
            this.isDragging = true;
            this.dragTarget = 'pendulum';
            this.canvas.style.cursor = 'grabbing';
        } else if (this.isPointInCart(x, y, cartX, cartY)) {
            this.isDragging = true;
            this.dragTarget = 'cart';
            this.canvas.style.cursor = 'grabbing';
        }

        this.lastMouseX = x;
        this.lastMouseY = y;
    }

    onMouseMove(e) {
        const pos = this.getMousePos(e);
        this.handleInteraction(pos.x, pos.y);
    }

    onTouchMove(e) {
        e.preventDefault();
        const pos = this.getTouchPos(e);
        this.handleInteraction(pos.x, pos.y);
    }

    handleInteraction(x, y) {
        if (!this.isDragging) {
            // Change le curseur selon l'element survole
            const cartX = this.centerX + this.currentX * this.scale;
            const cartY = this.centerY;

            const poleLength = this.currentL * this.scale;
            const poleEndX = cartX + poleLength * Math.sin(this.currentTheta);
            const poleEndY = cartY - poleLength * Math.cos(this.currentTheta);

            if (this.isPointInPendulum(x, y, poleEndX, poleEndY) ||
                this.isPointInCart(x, y, cartX, cartY)) {
                this.canvas.style.cursor = 'grab';
            } else {
                this.canvas.style.cursor = 'default';
            }
            return;
        }

        const deltaX = x - this.lastMouseX;
        const deltaY = y - this.lastMouseY;

        if (this.dragTarget === 'cart') {
            // Deplacement du chariot
            const deltaPos = deltaX / this.scale;
            if (this.onCartDrag) {
                this.onCartDrag(deltaPos);
            }
        } else if (this.dragTarget === 'pendulum') {
            // Rotation du pendule
            const cartX = this.centerX + this.currentX * this.scale;
            const cartY = this.centerY;

            // Calcule l'angle base sur la position de la souris
            const dx = x - cartX;
            const dy = cartY - y; // Inverse car y croit vers le bas
            const newTheta = Math.atan2(dx, dy);

            if (this.onPendulumDrag) {
                this.onPendulumDrag(newTheta);
            }
        }

        this.lastMouseX = x;
        this.lastMouseY = y;
    }

    onMouseUp(e) {
        this.endInteraction();
    }

    onTouchEnd(e) {
        e.preventDefault();
        this.endInteraction();
    }

    endInteraction() {
        this.isDragging = false;
        this.dragTarget = null;
        this.canvas.style.cursor = 'default';
    }

    draw(x, theta, l, targetX = null) {
        // Stocke les valeurs actuelles pour l'interactivite
        this.currentX = x;
        this.currentTheta = theta;
        this.currentL = l;

        // Efface le canvas
        this.ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);

        // Dessine le rail
        this.drawRail();

        // Position du chariot en pixels
        const cartX = this.centerX + x * this.scale;

        // Dessine le chariot
        this.drawCart(cartX, this.centerY);

        // Dessine le pendule
        this.drawPole(cartX, this.centerY, theta, l);

        // Dessine la cible si mode projet
        if (targetX !== null) {
            this.drawTarget(this.centerX + targetX * this.scale, this.centerY);
        }

        // Affiche les valeurs
        this.drawValues(x, theta);

        // Affiche les instructions d'interactivite
        this.drawInstructions();
    }

    drawRail() {
        this.ctx.strokeStyle = '#34495e';
        this.ctx.lineWidth = 3;
        this.ctx.beginPath();
        this.ctx.moveTo(50, this.railY);
        this.ctx.lineTo(this.canvas.width - 50, this.railY);
        this.ctx.stroke();
    }

    drawCart(x, y) {
        // Surbrillance si survole ou en cours de deplacement
        const isHighlighted = this.dragTarget === 'cart' ||
            (!this.isDragging && this.isPointInCart(this.lastMouseX, this.lastMouseY, x, y));

        this.ctx.fillStyle = isHighlighted ? '#2980b9' : '#3498db';
        this.ctx.fillRect(
            x - this.cartWidth / 2,
            y - this.cartHeight / 2,
            this.cartWidth,
            this.cartHeight
        );

        // Contour si survole
        if (isHighlighted) {
            this.ctx.strokeStyle = '#1abc9c';
            this.ctx.lineWidth = 2;
            this.ctx.strokeRect(
                x - this.cartWidth / 2,
                y - this.cartHeight / 2,
                this.cartWidth,
                this.cartHeight
            );
        }

        // Roues
        this.ctx.fillStyle = '#2c3e50';
        const wheelRadius = 8;

        // Roue gauche
        this.ctx.beginPath();
        this.ctx.arc(x - this.cartWidth / 3, y + this.cartHeight / 2, wheelRadius, 0, 2 * Math.PI);
        this.ctx.fill();

        // Roue droite
        this.ctx.beginPath();
        this.ctx.arc(x + this.cartWidth / 3, y + this.cartHeight / 2, wheelRadius, 0, 2 * Math.PI);
        this.ctx.fill();
    }

    drawPole(cartX, cartY, theta, l) {
        const poleLength = l * this.scale;
        const poleEndX = cartX + poleLength * Math.sin(theta);
        const poleEndY = cartY - poleLength * Math.cos(theta);

        // Tige du pendule
        this.ctx.strokeStyle = '#e74c3c';
        this.ctx.lineWidth = this.poleWidth;
        this.ctx.beginPath();
        this.ctx.moveTo(cartX, cartY);
        this.ctx.lineTo(poleEndX, poleEndY);
        this.ctx.stroke();

        // Surbrillance de la masse si survole ou en cours de deplacement
        const isHighlighted = this.dragTarget === 'pendulum' ||
            (!this.isDragging && this.isPointInPendulum(this.lastMouseX, this.lastMouseY, poleEndX, poleEndY));

        // Masse du pendule
        this.ctx.fillStyle = isHighlighted ? '#c0392b' : '#e74c3c';
        this.ctx.beginPath();
        this.ctx.arc(poleEndX, poleEndY, 12, 0, 2 * Math.PI);
        this.ctx.fill();

        // Contour si survole
        if (isHighlighted) {
            this.ctx.strokeStyle = '#1abc9c';
            this.ctx.lineWidth = 2;
            this.ctx.beginPath();
            this.ctx.arc(poleEndX, poleEndY, 12, 0, 2 * Math.PI);
            this.ctx.stroke();
        }

        // Axe de rotation
        this.ctx.fillStyle = '#2c3e50';
        this.ctx.beginPath();
        this.ctx.arc(cartX, cartY, 6, 0, 2 * Math.PI);
        this.ctx.fill();
    }

    drawTarget(x, y) {
        this.ctx.strokeStyle = '#27ae60';
        this.ctx.lineWidth = 2;
        this.ctx.setLineDash([5, 5]);

        // Ligne verticale de la cible
        this.ctx.beginPath();
        this.ctx.moveTo(x, y - 40);
        this.ctx.lineTo(x, y + 40);
        this.ctx.stroke();

        this.ctx.setLineDash([]);
    }

    drawValues(x, theta) {
        this.ctx.fillStyle = '#2c3e50';
        this.ctx.font = '14px monospace';

        const xText = `x = ${x.toFixed(3)} m`;
        const thetaText = `θ = ${theta.toFixed(3)} rad (${(theta * 180 / Math.PI).toFixed(1)}°)`;

        this.ctx.fillText(xText, 10, 25);
        this.ctx.fillText(thetaText, 10, 45);
    }

    drawInstructions() {
        this.ctx.fillStyle = '#7f8c8d';
        this.ctx.font = '12px sans-serif';
        this.ctx.textAlign = 'right';

        const instructions = [
            'Cliquez et glissez:',
            '• Chariot pour le deplacer',
            '• Masse pour changer l\'angle'
        ];

        instructions.forEach((text, index) => {
            this.ctx.fillText(text, this.canvas.width - 10, this.canvas.height - 40 + index * 15);
        });

        this.ctx.textAlign = 'left';
    }
}

// ===== MODULE PLOTS (inchange) =====
class Plotter {
    constructor(canvasId, title, unit, color, timeWindow = 15) {
        this.canvas = document.getElementById(canvasId);
        this.ctx = this.canvas.getContext('2d');
        this.title = title;
        this.unit = unit;
        this.color = color;
        this.timeWindow = timeWindow;

        this.data = [];
        this.maxDataPoints = timeWindow * 100;
    }

    addPoint(time, value) {
        this.data.push({ t: time, y: value });

        while (this.data.length > 0 && this.data[0].t < time - this.timeWindow) {
            this.data.shift();
        }
    }

    draw() {
        if (this.data.length === 0) return;

        this.ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);

        const margin = { top: 30, right: 20, bottom: 40, left: 60 };
        const plotWidth = this.canvas.width - margin.left - margin.right;
        const plotHeight = this.canvas.height - margin.top - margin.bottom;

        const timeRange = this.timeWindow;
        const values = this.data.map(d => d.y);
        const yMin = Math.min(...values) * 1.1;
        const yMax = Math.max(...values) * 1.1;
        const yRange = yMax - yMin || 1;

        this.drawAxes(margin, plotWidth, plotHeight, yMin, yMax);
        this.drawCurve(margin, plotWidth, plotHeight, timeRange, yMin, yRange);

        this.ctx.fillStyle = '#2c3e50';
        this.ctx.font = 'bold 16px sans-serif';
        this.ctx.textAlign = 'center';
        this.ctx.fillText(this.title, this.canvas.width / 2, 20);
    }

    drawAxes(margin, plotWidth, plotHeight, yMin, yMax) {
        this.ctx.strokeStyle = '#34495e';
        this.ctx.lineWidth = 1;

        this.ctx.beginPath();
        this.ctx.moveTo(margin.left, margin.top);
        this.ctx.lineTo(margin.left, margin.top + plotHeight);
        this.ctx.stroke();

        this.ctx.beginPath();
        this.ctx.moveTo(margin.left, margin.top + plotHeight);
        this.ctx.lineTo(margin.left + plotWidth, margin.top + plotHeight);
        this.ctx.stroke();

        this.ctx.fillStyle = '#7f8c8d';
        this.ctx.font = '12px sans-serif';
        this.ctx.textAlign = 'right';

        for (let i = 0; i <= 4; i++) {
            const y = margin.top + (plotHeight * i) / 4;
            const value = yMax - (yMax - yMin) * i / 4;
            this.ctx.fillText(value.toFixed(2), margin.left - 5, y + 4);
        }

        this.ctx.save();
        this.ctx.translate(15, margin.top + plotHeight / 2);
        this.ctx.rotate(-Math.PI / 2);
        this.ctx.textAlign = 'center';
        this.ctx.fillText(this.unit, 0, 0);
        this.ctx.restore();

        this.ctx.textAlign = 'center';
        const currentTime = this.data.length > 0 ? this.data[this.data.length - 1].t : 0;

        for (let i = 0; i <= 5; i++) {
            const x = margin.left + (plotWidth * i) / 5;
            const time = currentTime - this.timeWindow + (this.timeWindow * i) / 5;
            this.ctx.fillText(time.toFixed(1), x, margin.top + plotHeight + 20);
        }

        this.ctx.fillText('Temps (s)', margin.left + plotWidth / 2, this.canvas.height - 5);
    }

    drawCurve(margin, plotWidth, plotHeight, timeRange, yMin, yRange) {
        if (this.data.length < 2) return;

        this.ctx.strokeStyle = this.color;
        this.ctx.lineWidth = 2;
        this.ctx.beginPath();

        const currentTime = this.data[this.data.length - 1].t;

        this.data.forEach((point, index) => {
            const x = margin.left + plotWidth * (point.t - (currentTime - timeRange)) / timeRange;
            const y = margin.top + plotHeight - (plotHeight * (point.y - yMin) / yRange);

            if (index === 0) {
                this.ctx.moveTo(x, y);
            } else {
                this.ctx.lineTo(x, y);
            }
        });

        this.ctx.stroke();
    }

    clear() {
        this.data = [];
    }
}

// ===== MODULE PERFORMANCE (inchange) =====
class PerformanceMonitor {
    constructor() {
        this.reset();
    }

    reset() {
        this.startTime = null;
        this.stabilizationTime = null;
        this.riseTime = null;
        this.maxAngle = 0;
        this.steadyStateError = 0;
        this.mode = 'strict';

        this.angleHistory = [];
        this.positionHistory = [];
        this.timeHistory = [];
    }

    setMode(mode) {
        this.mode = mode;
    }

    update(time, x, theta, targetX = 0, targetTheta = 0) {
        if (this.startTime === null) {
            this.startTime = time;
        }

        const t = time - this.startTime;

        this.timeHistory.push(t);
        this.angleHistory.push(Math.abs(theta));
        this.positionHistory.push(x);

        while (this.timeHistory.length > 0 && this.timeHistory[0] < t - 10) {
            this.timeHistory.shift();
            this.angleHistory.shift();
            this.positionHistory.shift();
        }

        this.maxAngle = Math.max(this.maxAngle, Math.abs(theta));

        this.checkStabilization(t, theta, targetTheta);

        if (this.mode === 'project') {
            this.checkRiseTime(t, x, targetX);
            this.calculateSteadyStateError(x, targetX);
        }

        return this.getMetrics();
    }

    checkStabilization(t, theta, targetTheta) {
        const threshold = this.mode === 'strict' ? 0.05 : 0.35;
        const windowSize = 1.0;

        if (this.angleHistory.length < 2) return;

        const recentAngles = this.angleHistory.slice(-Math.floor(windowSize * 100));
        const isStable = recentAngles.every(angle => angle <= threshold);

        if (isStable && this.stabilizationTime === null && t > windowSize) {
            this.stabilizationTime = t;
        }
    }

    checkRiseTime(t, x, targetX) {
        if (this.riseTime !== null) return;

        const threshold = 0.9 * Math.abs(targetX);

        if (Math.abs(x) >= threshold) {
            this.riseTime = t;
        }
    }

    calculateSteadyStateError(x, targetX) {
        if (this.timeHistory.length < 300) return;

        const recentPositions = this.positionHistory.slice(-300);
        const avgPosition = recentPositions.reduce((a, b) => a + b, 0) / recentPositions.length;

        this.steadyStateError = Math.abs(avgPosition - targetX) / Math.abs(targetX) * 100;
    }

    getMetrics() {
        const criteria = this.mode === 'strict' ?
            { stabilityTime: 5, angleThreshold: 0.05 } :
            { stabilityTime: 5, angleThreshold: 0.35, riseTime: 0.5, errorThreshold: 2 };

        return {
            stabilityTime: this.stabilizationTime,
            stabilityStatus: this.getStatus(this.stabilizationTime, criteria.stabilityTime, true),

            maxAngle: this.maxAngle,
            angleStatus: this.getStatus(this.maxAngle, criteria.angleThreshold, false),

            riseTime: this.riseTime,
            riseStatus: this.mode === 'project' ?
                this.getStatus(this.riseTime, criteria.riseTime, true) : null,

            steadyError: this.steadyStateError,
            errorStatus: this.mode === 'project' ?
                this.getStatus(this.steadyStateError, criteria.errorThreshold, false) : null,

            systemStatus: this.getSystemStatus(),
            overallStatus: this.getOverallStatus()
        };
    }

    getStatus(value, threshold, lowerIsBetter) {
        if (value === null) return null;

        const isGood = lowerIsBetter ? value <= threshold : value >= threshold;
        return isGood ? 'success' : 'danger';
    }

    getSystemStatus() {
        if (this.maxAngle > Math.PI / 2) return 'Chute detectee';
        if (this.stabilizationTime !== null) return 'Stabilise';
        return 'En cours';
    }

    getOverallStatus() {
        if (this.maxAngle > Math.PI / 2) return 'danger';

        const metrics = this.getMetrics();
        if (metrics.stabilityStatus === 'success' &&
            metrics.angleStatus === 'success' &&
            (this.mode === 'strict' ||
             (metrics.riseStatus === 'success' && metrics.errorStatus === 'success'))) {
            return 'success';
        }

        return 'warning';
    }
}

// ===== APPLICATION PRINCIPALE AVEC INTERACTIVITE =====
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
        this.manualControl = false; // Mode controle manuel

        this.setupCallbacks();
        this.setupInteractivity();
        this.reset();

        // Demarre automatiquement la simulation pour montrer l'animation
        setTimeout(() => {
            this.play();
        }, 1000);
    }

    setupCallbacks() {
        this.ui.onPlay = () => this.play();
        this.ui.onPause = () => this.pause();
        this.ui.onReset = () => this.reset();
        this.ui.onImpulse = () => this.applyImpulse();
        this.ui.onAutoTune = () => this.autoTune();
        this.ui.onPresetChange = () => this.applyPreset();
        this.ui.onScenarioChange = () => this.loadScenario();
    }

    setupInteractivity() {
        // Callback pour le deplacement du chariot
        this.animation.onCartDrag = (deltaPos) => {
            const newX = this.physics.x + deltaPos;
            // Limite le deplacement
            const limitedX = Math.max(-2, Math.min(2, newX));
            this.physics.setState(limitedX, 0, this.physics.theta, this.physics.theta_dot);
            this.manualControl = true;

            // Remet le controle automatique apres un delai plus court
            setTimeout(() => {
                this.manualControl = false;
            }, 500);
        };

        // Callback pour la rotation du pendule
        this.animation.onPendulumDrag = (newTheta) => {
            // Limite l'angle
            const limitedTheta = Math.max(-Math.PI/2, Math.min(Math.PI/2, newTheta));
            this.physics.setState(this.physics.x, this.physics.x_dot, limitedTheta, 0);
            this.manualControl = true;

            // Remet le controle automatique apres un delai plus court
            setTimeout(() => {
                this.manualControl = false;
            }, 500);
        };
    }

    play() {
        if (!this.isRunning) {
            this.isRunning = true;
            this.lastTime = performance.now();
            this.gameLoop();
        }
    }

    pause() {
        this.isRunning = false;
        if (this.animationId) {
            cancelAnimationFrame(this.animationId);
        }
    }

    reset() {
        this.pause();

        const scenario = this.ui.scenarioSelect.value;
        let initialState;

        if (scenario === 'impulse') {
            // Position proche de l'equilibre pour tester l'impulsion
            initialState = [0, 0, 0.01, 0];
        } else if (scenario === 'position') {
            // Position pendante pour le mode projet
            initialState = [0, 0, 0.1, 0];
        } else {
            // Position pendante par defaut - animation spectaculaire
            initialState = [0, 0, 0.1, 0];
        }

        this.physics.reset(initialState);
        this.controller.reset();
        this.performance.reset();
        this.performance.setMode(this.ui.mode);
        this.manualControl = false;

        // Reset des integrales du controleur
        this.thetaIntegral = 0;
        this.xIntegral = 0;

        this.anglePlot.clear();
        this.positionPlot.clear();

        this.updateSimulation();

        // Demarre automatiquement la simulation pour l'animation
        if (scenario !== 'impulse') {
            setTimeout(() => {
                this.play();
            }, 500);
        }
    }

    applyImpulse() {
        this.physics.applyImpulse(1.0);
    }

    autoTune() {
        const gains = this.controller.autoTune();
        this.ui.setPIDParams(gains);
    }

    applyPreset() {
        const preset = this.ui.pidPreset.value;
        if (preset !== 'custom') {
            const gains = PIDController.getPreset(preset);
            this.ui.setPIDParams(gains);
        }
    }

    loadScenario() {
        this.reset();

        const scenario = this.ui.scenarioSelect.value;
        if (scenario === 'impulse') {
            // Pour l'impulsion, on attend que l'utilisateur demarre
            setTimeout(() => {
                setTimeout(() => this.applyImpulse(), 500);
            }, 100);
        } else if (scenario === 'position') {
            this.ui.setMode('project');
            // La simulation demarre automatiquement via reset()
        }
        // Pour 'nominal', la simulation demarre automatiquement via reset()
    }

    gameLoop() {
        if (!this.isRunning) return;

        const currentTime = performance.now();
        const deltaTime = (currentTime - this.lastTime) / 1000;

        if (deltaTime >= this.physics.dt) {
            this.update();
            this.lastTime = currentTime;
        }

        this.render();
        this.animationId = requestAnimationFrame(() => this.gameLoop());
    }

    update() {
        this.updateParameters();

        // Applique toujours le controle, meme en mode manuel (avec force reduite)
        let force = 0;
        if (!this.manualControl) {
            force = this.calculateControl();
        } else {
            // Mode manuel : controle reduit mais toujours actif
            force = this.calculateControl() * 0.3; // Force reduite en mode manuel
        }

        const simParams = this.ui.getSimulationParams();
        this.physics.dt = simParams.dt;
        this.physics.step(force, simParams.integrator);

        // Detection de problemes mais sans reset automatique
        // Le systeme doit continuer a essayer de se stabiliser
        // if (this.physics.hasFallen()) {
        //     console.log('Systeme instable detecte - le controle continue');
        // }

        const targetX = this.ui.mode === 'project' ? 0.2 : 0;
        const metrics = this.performance.update(
            this.physics.t,
            this.physics.x,
            this.physics.theta,
            targetX,
            0
        );

        this.anglePlot.addPoint(this.physics.t, this.physics.theta);
        this.positionPlot.addPoint(this.physics.t, this.physics.x);

        this.ui.updateIndicators(metrics);
    }

    updateParameters() {
        const pidParams = this.ui.getPIDParams();
        this.controller.setGains(pidParams.kp, pidParams.ki, pidParams.kd);

        const sysParams = this.ui.getSystemParams();
        this.physics.updateParams(sysParams);
    }

    calculateControl() {
        const theta = this.physics.theta;
        const theta_dot = this.physics.theta_dot;
        const x = this.physics.x;
        const x_dot = this.physics.x_dot;

        // Normalise l'angle entre -PI et PI (comme dans le modele de reference)
        let normalizedTheta = theta;
        while (normalizedTheta > Math.PI) normalizedTheta -= 2 * Math.PI;
        while (normalizedTheta < -Math.PI) normalizedTheta += 2 * Math.PI;

        if (this.ui.mode === 'strict') {
            // Controleur PD simple et efficace (formule du modele de reference)
            // F = -Kp*θ - Kd*θ˙
            // Efficace pres de l'equilibre

            const pidGains = this.ui.getPIDParams();
            const Kp = pidGains.kp;
            const Kd = pidGains.kd;

            // Force PD principale
            const F_pd = -Kp * normalizedTheta - Kd * theta_dot;

            // Petit terme de rappel de position pour eviter la derive
            const F_position = -1.0 * x - 1.5 * x_dot;

            return F_pd + F_position;

        } else {
            // Mode projet - suivi de position avec PD
            const pidGains = this.ui.getPIDParams();

            // Force pour l'angle
            const F_angle = -pidGains.kp * normalizedTheta - pidGains.kd * theta_dot;

            // Force pour la position (suivi de 0.2m)
            const positionError = 0.2 - x;
            const F_position = 10 * positionError - 3 * x_dot;

            return F_angle + F_position;
        }
    }

    render() {
        const targetX = this.ui.mode === 'project' ? 0.2 : null;
        this.animation.draw(
            this.physics.x,
            this.physics.theta,
            this.physics.l,
            targetX
        );

        this.anglePlot.draw();
        this.positionPlot.draw();
    }

    updateSimulation() {
        this.render();
    }
}

// ===== INITIALISATION =====
document.addEventListener('DOMContentLoaded', () => {
    const app = new InvertedPendulumApp();
    window.pendulumApp = app;
});