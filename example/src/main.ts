import { stateToKeplerian, type Keplerian, type State } from "sugarfree-spearmint";

import {
  AxesHelper,
  BufferGeometry,
  EllipseCurve,
  Group,
  IcosahedronGeometry,
  LineBasicMaterial,
  LineLoop,
  Mesh,
  MeshBasicMaterial,
  PerspectiveCamera,
  Scene,
  SphereGeometry,
  Vector3,
  WebGLRenderer,
} from "three";
import { OrbitControls } from "three/examples/jsm/Addons.js";

const renderer = new WebGLRenderer();
const scene = new Scene();
const camera = new PerspectiveCamera();
const controls = new OrbitControls(camera, renderer.domElement);
controls.autoRotate = true;
controls.autoRotateSpeed = -10;

const scale = 0.001;
const earthRadius = 6371 * scale;
const earthMu = 398600.44 * scale;

const issState: State = {
  position: new Vector3(-2538.5121656234, -4211.67676741827, 4680.2482585403).multiplyScalar(scale).toArray(),
  velocity: new Vector3(4.19261976686944, -5.73481146139295, -2.88408692090738).toArray(),
};

const issKeplerian = stateToKeplerian(earthMu, issState);

const earthMesh = new Mesh(
  new SphereGeometry(earthRadius),
  new MeshBasicMaterial({ color: "#1A428A", wireframe: true }),
);
const issMesh = new Mesh(new IcosahedronGeometry(0.2, 1), new MeshBasicMaterial({ color: "red" }));

const issOrbitMesh = createOrbitMesh(issKeplerian);

function start() {
  scene.add(earthMesh);
  issMesh.position.set(...issState.position);
  scene.add(issMesh);
  scene.add(issOrbitMesh);

  scene.add(new AxesHelper(10));
  // camera.position.set(-20, 20, 20);
  camera.position.setZ(30).setY(10);

  resize();
  window.addEventListener("resize", resize);
  document.querySelector("#app")?.appendChild(renderer.domElement);
  renderer.setAnimationLoop(callback);
}

function callback() {
  controls.update();
  renderer.render(scene, camera);
}

function resize() {
  renderer.setSize(window.innerWidth, window.innerHeight);
  camera.aspect = window.innerWidth / window.innerHeight;
  camera.updateProjectionMatrix();
  controls.update();
}

export function createOrbitMesh(keplerian: Keplerian, color = 0xffffff) {
  const { semiMajorAxis: a, eccentricity: e, inclination: i, raan, argumentOfPeriapsis: aop } = keplerian;

  // 1. Calculate ellipse dimensions
  const b = a * Math.sqrt(1 - e ** 2); // Semi-minor axis
  const focusOffset = a * e; // Distance from center to focus

  // 2. Create the 2D Ellipse Curve in the XY plane
  // Note: Focus is at (0,0), so we shift the center to (-focusOffset, 0)
  const curve = new EllipseCurve(
    -focusOffset,
    0, // Center offset
    a,
    b, // xRadius, yRadius
    0,
    2 * Math.PI, // Start, End angle
    false, // Clockwise
    0, // Rotation
  );

  const points = curve.getPoints(128);
  const geometry = new BufferGeometry().setFromPoints(points);
  const material = new LineBasicMaterial({ color });
  const orbitLine = new LineLoop(geometry, material);

  // 3. Apply 3D Orientation using Nested Wrappers
  // Order: RAAN (Z) -> Inclination (X) -> Argument of Periapsis (Z)
  const raanGroup = new Group();
  const inclinationGroup = new Group();
  const aopGroup = new Group();

  raanGroup.add(inclinationGroup);
  inclinationGroup.add(aopGroup);
  aopGroup.add(orbitLine);

  // Apply rotations (Three.js uses Radians)
  raanGroup.rotation.z = raan;
  inclinationGroup.rotation.x = i;
  aopGroup.rotation.z = aop;

  return raanGroup;
}

start();
