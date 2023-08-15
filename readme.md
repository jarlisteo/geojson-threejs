# Libreria que convierte un archivo GeoJSON a una visualizacion de ThreeJS

## Comandos de la libreria:

- npm install -> instala dependencias
- npm vite -> corre vite en modo dev para previsualizar la libreria
- npm vite build -> hace un build de la libreria en /dist

## Build

El build es una libreria javascript que contiene todo lo necesario para montarlo en cualquier lugar

El uso es muy sencillo se muestra a continuacion

```javascript
<script type="module">import LatbitGeojson from "./main.js"; new LatbitGeojson(0);</script>
```

Esto crea en el body un visualizador threejs para mostrar el GEOJSON del id, en el ejemplo se usa ID = 0 que es el geojson de prueba.
