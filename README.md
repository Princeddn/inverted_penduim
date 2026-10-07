# Pendule inversé sur chariot

Simulation web temps réel d'un pendule inversé stabilisé par retour d'état.

## Lancer

Ouvrir `index.html` dans un navigateur (aucune dépendance, aucun serveur).

## Contenu

- `index.html` : interface (contrôles, indicateurs, scène, courbes).
- `script.js` : physique (Lagrange + RK4/Euler), contrôleur, auto-réglage par placement de pôles, animation, tracés, indicateurs.
- `style.css` : mise en forme.

## Utilisation

- **Départ en bas** : le pendule pend au repos ; « Lancer » le redresse (pompage d’énergie) puis le stabilise.
- **Stabilisation** : ramène θ → 0 et x → 0.
- **Suivi x = 0.2 m** : déplace le chariot vers 0.2 m en gardant le pendule debout.
- **Impulsion** : 1 N·s appliqué au chariot.
- **Auto-réglage** : recalcule les gains pour les paramètres physiques actuels.
- Glisser le chariot ou la masse à la souris ; flèches ← / → pour pousser le chariot.

## Modèle physique (banc de laboratoire)

- Chariot de 1 kg sur un rail de ±0,6 m terminé par des butées en caoutchouc (le choc est transmis au pendule).
- Pendule : tige homogène de 0,6 m et 250 g.
- Moteur à courant continu : F = αV − βẋ (α = 1,72 N/V, β = 7,7 N·s/m, |V| ≤ 10 V), soit 17 N max à l'arrêt et environ 2,1 m/s max.
- Calculateur échantillonné (5 ms), encodeurs quantifiés et bruités, vitesses estimées, frottement sec.
- La souris agit comme une main (lien élastique) sur le chariot ou le bout de la tige ; la simulation ne s'arrête jamais.
