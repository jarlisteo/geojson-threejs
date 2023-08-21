import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";

import exampleGeojson from "./example1.json" assert { type: "json" };

var scene, camera, renderer, controls;
var selectedId = 0;

const raycaster = new THREE.Raycaster();
const pointer = new THREE.Vector2();
const colors = ["#00a5e3", "#8dd7bf", "#6c88c4", "#ffa23a", "#ffd872", "#ff3d18", "#ff6446", "#ff8b74"];

export default class LatbitGeojson {
  constructor(id, endpoint = "http://108.181.190.199/test/get_unidades.php?id=") {
    this.id = id;
    this.geojson = exampleGeojson;
    this.endpoint = endpoint;

    if (id != 0) {
      this.geojson = obtenerGeoJSON(endpoint, id).then((json) => {
        this.geojson = json;
        this.render();
      });
      return;
    }

    return this.render();
  }

  render() {
    renderer = new THREE.WebGLRenderer({ antialias: true });
    renderer.shadowMap.enabled = true;
    renderer.setSize(window.innerWidth, window.innerHeight);

    document.body.appendChild(renderer.domElement);

    scene = new THREE.Scene();
    camera = new THREE.PerspectiveCamera(45, window.innerWidth / window.innerHeight, 0.1, 1000);

    controls = new OrbitControls(camera, renderer.domElement);

    this.x_values = [];
    this.z_values = [];

    this.axis1_min = 0;
    this.axis3_min = 0;
    this.axis1_max = 0;
    this.axis3_max = 0;

    var scale = 5000;
    var height_scale = 0.05;

    this.calculateMinAndMax();

    renderer.setClearColor(0xffffff, 0);

    const light = new THREE.DirectionalLight("white", 2);
    light.castShadow = true;
    light.shadow.mapSize.width = 10240;
    light.shadow.mapSize.height = 10240;
    light.shadow.radius = 10;

    const alight = new THREE.AmbientLight(0x404040, 20); // soft white light
    scene.add(alight);

    light.position.set(1, 2, 1);
    scene.add(light);

    this.geojson.features.forEach((feature) => {
      feature.geometry.coordinates.forEach((polygon) => {
        this.convertToPlaneCoords(polygon[0]);
      });
      let nivel = feature.properties.nivel;
      let height = feature.properties.altura * height_scale;
      let vertices = this.getVertex(this.x_values, height, this.z_values);

      let planePoints = vertices.map((v) => {
        return new THREE.Vector2(v.x, v.z);
      });

      const extrudeSettings = {
        steps: 1,
        depth: height,
        bevelEnabled: false,
        bevelThickness: 1,
        bevelSize: 1,
        bevelOffset: 0,
        bevelSegments: 1,
      };

      const shape = new THREE.Shape(planePoints);
      const geometry = new THREE.ExtrudeGeometry(shape, extrudeSettings);
      const shapeMesh = new THREE.Mesh(
        geometry,
        new THREE.MeshPhysicalMaterial({ color: colors[feature.properties.nivel], side: THREE.DoubleSide, transparent: true, opacity: 0.8 })
      );

      shapeMesh.scale.set(scale, scale, 1);
      shapeMesh.rotateX(Math.PI / 2);
      shapeMesh.position.set(0, height * nivel, 0);
      shapeMesh.castShadow = true;
      shapeMesh.id_piso = feature.properties.id_piso;
      shapeMesh.nivel = nivel;
      scene.add(shapeMesh);

      this.clearArrays();
    });

    const planeGeometry = new THREE.PlaneGeometry(100, 100);
    const planeMaterial = new THREE.MeshStandardMaterial({ color: "white", side: THREE.DoubleSide });
    const plane = new THREE.Mesh(planeGeometry, planeMaterial);
    plane.receiveShadow = true;
    plane.rotateX(Math.PI / 2);
    scene.add(plane);

    camera.position.set(0.5, 1, -1.5);

    controls.update();
    animate();
  }

  convertToPlaneCoords(vertexs) {
    vertexs.forEach((vertex) => {
      let lon = vertex[0];
      let lat = vertex[1];
      this.x_values.push(lat);
      this.z_values.push(lon);
    });
  }

  calculateMinAndMax() {
    let allX = [];
    let allZ = [];
    let allNivels = [];

    this.geojson.features.forEach((feature) => {
      feature.properties.nivel;
      allNivels.push(feature.properties.nivel);
      feature.geometry.coordinates.forEach((polygon) => {
        polygon[0].forEach((vertex) => {
          allX.push(vertex[1]);
          allZ.push(vertex[0]);
        });
      });
    });

    this.axis1_min = Math.min(...allX);
    this.axis3_min = Math.min(...allZ);
    this.axis1_max = Math.max(...allX);
    this.axis3_max = Math.max(...allZ);
  }

  getVertex(values_axis1, values_axis2, values_axis3) {
    let points = [];
    for (var i = 0; i < values_axis1.length; i++) {
      let x = values_axis1[i] - (this.axis1_min + this.axis1_max) / 2;
      let z = values_axis3[i] - (this.axis3_min + this.axis3_max) / 2;
      points.push(new THREE.Vector3(x, values_axis2, z));
    }
    return points;
  }

  clearArrays() {
    this.x_values = [];
    this.z_values = [];
  }
}

async function obtenerGeoJSON(endpoint, ID) {
  try {
    let res = await fetch(endpoint + ID);
    let json = await res.json();
    return json;
  } catch (error) {
    console.log("Error al obtener el GeoJSON desde el endpoint, verifique que el endpoint sea correcto");
  }
}

function animate() {
  requestAnimationFrame(animate);
  controls.update();
  render();
}

function render() {
  handleSelected();
  handleSelectedMaterial();
  renderer.render(scene, camera);
}

function handleSelectedMaterial() {
  scene.children.forEach((child) => {
    if (child.nivel != undefined) {
      if (child.id == selectedId) {
        child.material = new THREE.MeshPhysicalMaterial({ color: "red", side: THREE.DoubleSide, transparent: true, opacity: 0.8 });
      } else {
        child.material = new THREE.MeshPhysicalMaterial({ color: colors[child.nivel], side: THREE.DoubleSide, transparent: true, opacity: 0.8 });
      }
    }
  });
}

function handleSelected() {
  if (pointer.x != 0 && pointer.y != 0) {
    raycaster.setFromCamera(pointer, camera);
    let intersects = raycaster.intersectObjects(scene.children);
    if (intersects.length > 0) {
      if (intersects[0].object.id_piso != undefined) {
        console.log("Nivel: " + intersects[0].object.nivel + " ID: " + intersects[0].object.id_piso);
        selectedId = intersects[0].object.id;
      }
    }
  }
}

function onPointerMove(event) {
  pointer.x = (event.clientX / window.innerWidth) * 2 - 1;
  pointer.y = -(event.clientY / window.innerHeight) * 2 + 1;
}

function resetPointerMove() {
  setTimeout(() => {
    pointer.x = 0;
    pointer.y = 0;
  }, 10);
}

window.addEventListener("mousedown", onPointerMove);
window.addEventListener("mouseup", resetPointerMove);
