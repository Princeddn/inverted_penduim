# Pendule inversé sur chariot

Simulation web temps réel d'un pendule inversé stabilisé par retour d'état.

## Lancer

Ouvrir `index.html` dans un navigateur (aucune dépendance, aucun serveur).

## Contenu

- `index.html` : interface (contrôles, indicateurs, scène, courbes).
- `script.js` : physique (Lagrange + RK4/Euler), contrôleur, auto-réglage par placement de pôles, animation, tracés, indicateurs.
- `style.css` : mise en forme.

## Utilisation

- **Stabilisation** : ramène θ → 0 et x → 0.
- **Suivi x = 0.2 m** : déplace le chariot vers 0.2 m en gardant le pendule debout.
- **Impulsion** : 1 N·s appliqué au chariot.
- **Auto-réglage** : recalcule les gains pour les paramètres physiques actuels.
- Glisser le chariot ou la masse à la souris ; flèches ← / → pour pousser le chariot.
