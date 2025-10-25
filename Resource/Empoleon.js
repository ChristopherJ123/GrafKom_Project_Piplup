import { Geometry } from "./Geometries.js";
import { LIBS } from "./Libs.js";
import { ModelNode } from "./ModelNode.js";

export class Empoleon {
    constructor(gl, renderer) {
        this.gl = gl;
        this.renderer = renderer;

        // NEW: Root of the scene graph
        this.rootNode = new ModelNode(gl);

        // This will be its static offset matrix (like Y+2.5)
        this.modelMatrix = LIBS.get_I4();

        // NEW: For future animations
        this.animatedNodes = {
            breathingNode: null, // <-- ADDED: Node to control overall Y position + animation
            leftHand: null,      // <-- ADDED: Node for left hand geometry
            rightHand: null      // <-- ADDED: Node for right hand geometry
        };
        this.baseTransforms = {
            leftHand: LIBS.get_I4(), // <-- ADDED: Base pose for left hand
            rightHand: LIBS.get_I4() // <-- ADDED: Base pose for right hand
        };

        this.initParts();
    }

    initParts() {
        const gl = this.gl;
        // Empoleon Colors
        const C = {
            BODY: [0.52, 0.80, 1.00], HEAD: [0.294, 0.541, 0.796], BEAK: [1.00, 0.84, 0.00], CAPE: [0.24, 0.42, 0.96],
            EYE_W: [1.00, 1.00, 1.00], BLACK: [0.00, 0.00, 0.00], FEET: [1.00, 0.65, 0.00], WHITE: [1.00, 1.00, 1.00],
            TAIL: [0.30, 0.54, 0.80], EMPO_BASE: [0.2, 0.247, 0.278], EMPO_LOWER_BODY: [0.18, 0.224, 0.247]
        };

        // Helper function to create a translation matrix using your libs.js functions
        const createTransform = (x, y, z) => {
            const m = LIBS.get_I4();
            LIBS.translateX(m, x);
            LIBS.translateY(m, y);
            LIBS.translateZ(m, z);
            return m;
        };

        const headTexture = this.renderer.loadTexture("Resource/empoleon_texture.png");

        // --- NEW HIERARCHY SETUP ---
        // 1. Create the invisible "breathing" node.
        const breathingNode = new ModelNode(gl);
        this.rootNode.addChild(breathingNode);
        this.animatedNodes.breathingNode = breathingNode;
        // --- END NEW HIERARCHY ---

        const body_profile = [
            // Duplicates removed for correct normals
            [0.0, 1.95, 0],
            [0.05, 1.95, 0],
            [0.25, 1.9, 0],
            [0.375, 1.8, 0],
            [0.425, 1.7, 0],
            [0.45, 1.6, 0],
            [0.475, 1.5, 0],
            [0.4875, 1.45, 0], // eyes stretching
            [0.5, 1.4, 0],
            [0.5125, 1.35, 0],
            [0.525, 1.3, 0],
            [0.5375, 1.25, 0],
            [0.55, 1.2, 0],   // Top point (neck)
            [0.56, 1.15, 0],
            [0.57, 1.1, 0],
            [0.59, 1.0, 0],
            [0.62, 0.9, 0],
            [0.65, 0.8, 0],
            [0.7, 0.7, 0],
            [0.75, 0.6, 0],
            [0.8, 0.5, 0],
            [0.85, 0.4, 0],
            [0.9, 0.3, 0],
            [0.95, 0.2, 0],
            [1, 0.1, 0],
            [1.05, 0, 0],
            [1.1, -0.1, 0],
            [1.15, -0.2, 0],
            [1.2, -0.3, 0], // Widest part of the belly
            [1.1875, -0.4, 0],
            [1.17, -0.5, 0],
            [1.14, -0.6, 0],
            [1.1, -0.7, 0],
            [1.05, -0.8, 0],
            [1.0, -0.9, 0],
            [0.9, -1.0, 0],
            [0.8, -1.1, 0],
            [0.55, -1.2, 0],
            [0.0, -1.3, 0]   // Bottom point
        ];

        // Define parts and their local transformations
        // --- NOTE: All transforms defined here should be RELATIVE to their parent ---
        // --- Since all parts are children of breathingNode, these are relative to (0,0,0) for now ---
        const partDefinitions = [
            // NEW Prinplup Body
            {
                geom: Geometry.generateLathe(body_profile, 30, C.BODY),
                trans: (() => {
                    const m = createTransform(0, 0.2, 0);
                    LIBS.rotateY(m, Math.PI / 2);
                    LIBS.scale(m, 1.2);
                    return m
                })(), // Slightly raise the body
                texture: headTexture
            },

            {
                geom: Geometry.generateTaperedShapeFromSpline(
                    [
                        [0.0, -0.6, -0.8], [0.0, -0.6, -1.2], [0.0, -0.6, -2.0]
                    ], 50, [0.2, 0.9], [0.01, 0.01], 20, C.HEAD
                ),
                trans: createTransform(0, 0.6, 0.1),
            },

            {
                geom: Geometry.generateTaperedShapeFromSpline(
                    [
                        [0.0, 0.2, 1.0], [0.0, 0.2, 1.2], [0.0, 0.2, 1.4]
                    ], 50, [0.15, 0.9], [0.01, 0.01], 20, C.HEAD
                ),
                trans: (() => {
                    const m = createTransform(0, -0.35, 0.1);
                    LIBS.rotateX(m, -Math.PI / 8);
                    return m
                })(),
            },

            // BEAK
            { geom: Geometry.generateTriBeak(0.7, 0.1, 0.5, C.BEAK), trans: (() => {
                    let m = createTransform(0, 1.79, 0.30);
                    LIBS.rotateZ(m, LIBS.degToRad(0));
                    LIBS.rotateY(m, LIBS.degToRad(180));
                    LIBS.rotateX(m, LIBS.degToRad(0));
                    return m;
                })()},

            { geom: Geometry.generateTriBeak(0.6, 0.1, 0.5, C.BEAK), trans: (() => {
                    let m = createTransform(0, 1.69, 0.30);
                    LIBS.rotateZ(m, LIBS.degToRad(180));
                    LIBS.rotateY(m, LIBS.degToRad(180));
                    LIBS.rotateX(m, LIBS.degToRad(0));
                    return m;
                })()},

            // CROWN
            { geom: Geometry.generateTriangularPrism(0.1, 0.1, 1.15, C.BEAK), trans: (() => {
                    let m = createTransform(0, 2.37, 0.45);
                    LIBS.rotateZ(m, LIBS.degToRad(180));
                    LIBS.rotateY(m, LIBS.degToRad(0));
                    LIBS.rotateX(m, LIBS.degToRad(-100));
                    return m;
                })()},

            { geom: Geometry.generateTriangularPrism(0.1, 0.1, 0.7, C.BEAK), trans: (() => {
                    let m = createTransform(0.43, 2.4, 0.20);
                    LIBS.rotateZ(m, LIBS.degToRad(210));
                    LIBS.rotateY(m, LIBS.degToRad(0));
                    LIBS.rotateX(m, LIBS.degToRad(-100));
                    return m;
                })()},

            { geom: Geometry.generateTriangularPrism(0.1, 0.1, 0.7, C.BEAK), trans: (() => {
                    let m = createTransform(-0.43, 2.4, 0.20);
                    LIBS.rotateZ(m, LIBS.degToRad(-210));
                    LIBS.rotateY(m, LIBS.degToRad(0));
                    LIBS.rotateX(m, LIBS.degToRad(-100));
                    return m;
                })()},

            { geom: Geometry.generateTriangularPrism(0.1, 0.1, 0.5, C.BEAK), trans: (() => {
                    let m = createTransform(0.30, 1.9, 0.37);
                    LIBS.rotateZ(m, LIBS.degToRad(-15));
                    LIBS.rotateY(m, LIBS.degToRad(33));
                    LIBS.rotateX(m, LIBS.degToRad(-125));
                    return m;
                })()},

            { geom: Geometry.generateTriangularPrism(0.1, 0.1, 0.5, C.BEAK), trans: (() => {
                    let m = createTransform(-0.30, 1.9, 0.37);
                    LIBS.rotateZ(m, LIBS.degToRad(15));
                    LIBS.rotateY(m, LIBS.degToRad(-33));
                    LIBS.rotateX(m, LIBS.degToRad(-125));
                    return m;
                })()},

            { geom: Geometry.generateComplexCone(0.1, 0, 0.2, 3, 20, C.BEAK), trans: (() => {
                    let m = createTransform(0, 3.045, 0.317);
                    LIBS.rotateZ(m, LIBS.degToRad(0));
                    LIBS.rotateY(m, LIBS.degToRad(30));
                    LIBS.rotateX(m, LIBS.degToRad(-10));
                    return m;
                })()},

            { geom: Geometry.generateComplexCone(0.1, 0, 0.2, 3, 20, C.BEAK), trans: (() => {
                    let m = createTransform(0, 2.85, 0.35);
                    LIBS.rotateZ(m, LIBS.degToRad(0));
                    LIBS.rotateY(m, LIBS.degToRad(210));
                    LIBS.rotateX(m, LIBS.degToRad(-190));
                    return m;
                })()},

            { geom: Geometry.generateComplexCone(0.1, 0, 0.2, 3, 20, C.BEAK), trans: (() => {
                    let m = createTransform(-0.422, 2.92, 0.087);
                    LIBS.rotateZ(m, LIBS.degToRad(0));
                    LIBS.rotateY(m, LIBS.degToRad(0));
                    LIBS.rotateX(m, LIBS.degToRad(-10));
                    return m;
                })()},

            { geom: Geometry.generateComplexCone(0.1, 0, 0.2, 3, 20, C.BEAK), trans: (() => {
                    let m = createTransform(-0.422, 2.73, 0.123);
                    LIBS.rotateZ(m, LIBS.degToRad(0));
                    LIBS.rotateY(m, LIBS.degToRad(0));
                    LIBS.rotateX(m, LIBS.degToRad(-190));
                    return m;
                })()},

            { geom: Geometry.generateComplexCone(0.1, 0, 0.2, 3, 20, C.BEAK), trans: (() => {
                    let m = createTransform(0.422, 2.92, 0.087);
                    LIBS.rotateZ(m, LIBS.degToRad(0));
                    LIBS.rotateY(m, LIBS.degToRad(60));
                    LIBS.rotateX(m, LIBS.degToRad(-10));
                    return m;
                })()},

            { geom: Geometry.generateComplexCone(0.1, 0, 0.2, 3, 20, C.BEAK), trans: (() => {
                    let m = createTransform(0.422, 2.73, 0.123);
                    LIBS.rotateZ(m, LIBS.degToRad(0));
                    LIBS.rotateY(m, LIBS.degToRad(60));
                    LIBS.rotateX(m, LIBS.degToRad(-190));
                    return m;
                })()},

            // KERAH Kiri
            {
                geom: Geometry.generateTaperedShapeFromSpline(
                    [
                        [0.0, 0.2, 1.4], [0.0, 0.2, 1.2], [0.0, 0.2, 1.0]
                    ], 50, [0.02, 0.5], [0.01, 0.01], 20, C.HEAD
                ),
                trans: (() => {
                    const m = createTransform(0.17, 1.72 , -0.85);
                    LIBS.rotateX(m, LIBS.degToRad(-15));
                    LIBS.rotateY(m, LIBS.degToRad(-15));
                    LIBS.rotateZ(m, LIBS.degToRad(40));
                    return m
                })(),
            },
            {
                geom: Geometry.generateTaperedShapeFromSpline(
                    [
                        [0.0, 0.2, 1.4], [0.0, 0.2, 1.2], [0.0, 0.2, 1.0]
                    ], 50, [0.02, 0.5], [0.01, 0.01], 20, C.HEAD
                ),
                trans: (() => {
                    const m = createTransform(-0.17, 1.72 , -0.85);
                    LIBS.rotateX(m, LIBS.degToRad(-15));
                    LIBS.rotateY(m, LIBS.degToRad(15));
                    LIBS.rotateZ(m, LIBS.degToRad(-40));
                    return m
                })(),
            },

            { geom: Geometry.generateSphere(0.1, 0.3, 0.1, 10, 10, C.BEAK), trans: (() => {
                    const m = createTransform(0.4, 1.25, -0.5);
                    LIBS.rotateX(m, Math.PI / 8);
                    return m
                })()},

            { geom: Geometry.generateSphere(0.1, 0.3, 0.1, 10, 10, C.BEAK), trans: (() => {
                    const m = createTransform(-0.4, 1.25, -0.5);
                    LIBS.rotateX(m, Math.PI / 8);
                    return m
                })()},


            // upper front body decor
            { geom: Geometry.generateSphere(0.08, 0.8, 0.08, 10, 10, C.HEAD), trans: (() => {
                    const m = createTransform(0, 0.8, 0.85);
                    LIBS.rotateX(m, -Math.PI / 8);
                    return m
                })()},

            { geom: Geometry.generateSphere(0.08, 0.8, 0.08, 10, 10, C.HEAD), trans: (() => {
                    const m = createTransform(0, 0.6, 0.95);
                    LIBS.rotateX(m, -Math.PI / 8);
                    return m
                })()},

            { geom: Geometry.generateSphere(0.1, 0.2, 0.15, 10, 10, C.HEAD), trans: (() => {
                    const m = createTransform(0, -0.5, 1.1);
                    LIBS.rotateX(m, 0.2);
                    return m
                })()},

            { geom: Geometry.generateSphere(0.1, 0.2, 0.15, 10, 10, C.HEAD), trans: (() => {
                    const m = createTransform(0, -0.6, 1.09);
                    LIBS.rotateX(m, 1.2);
                    return m
                })()},

            // --- Hands ---
            // MODIFIED: Added names and parentName (parent is breathingNode)
            {
                name: 'leftHand', // <-- ADDED
                parentName: 'breathingNode', // <-- ADDED
                geom: Geometry.generateSphere(0.2, 1.5, 0.5, 15, 15, C.EMPO_BASE),
                trans: (() => {
                    let m = createTransform(-1.2, 0.25, 0.2);
                    LIBS.rotateZ(m, LIBS.degToRad(-30));
                    LIBS.rotateX(m, LIBS.degToRad(-10));
                    return m;
                })()
            },
            {
                name: 'rightHand', // <-- ADDED
                parentName: 'breathingNode', // <-- ADDED
                geom: Geometry.generateSphere(0.2, 1.5, 0.5, 15, 15, C.EMPO_BASE),
                trans: (() => {
                    let m = createTransform(1.2, 0.25, 0.2);
                    LIBS.rotateZ(m, LIBS.degToRad(30));
                    LIBS.rotateX(m, LIBS.degToRad(-10));
                    return m;
                })()
            },

            // --- Outer Hands / Spikes ---
            // MODIFIED: Added parentName (parent is the corresponding hand)
            // Primary blue hand
            {
                parentName: 'leftHand', // <-- ADDED
                geom: Geometry.generateHalfEllipsoid(0.15, 1.5, 0.6, 15, 15, C.HEAD, 'lower'),
                trans: (() => {
                    // This transform MUST become relative to the parent hand later
                    let m = LIBS.get_I4();
                    // let m = createTransform(-1.225, 0.20, 0.2);
                    // LIBS.rotateZ(m, LIBS.degToRad(-30));
                    // LIBS.rotateX(m, LIBS.degToRad(-10));
                    return m;
                })()
            },
            {
                parentName: 'leftHand', // <-- ADDED
                geom: Geometry.generateComplexCone(0.2, 0, 0.4, 4, 20, C.HEAD),
                trans: (() => {
                    // This transform MUST become relative to the parent hand later
                    // let m = LIBS.get_I4();

                    let m = createTransform(0, 0, 0.55);
                    LIBS.rotateZ(m, LIBS.degToRad(180));
                    LIBS.rotateY(m, LIBS.degToRad(30));
                    LIBS.rotateX(m, LIBS.degToRad(-100));
                    return m;
                })()
            },
            // ... (Repeat adding parentName: 'leftHand' or 'rightHand' for all spikes and fingers) ...
            {
                parentName: 'leftHand', geom: Geometry.generateComplexCone(0.2, 0, 0.4, 4, 20, C.HEAD),
                trans: (() => {
                    // let m = LIBS.get_I4();
                    let m = createTransform(0, 0, -0.55);
                    LIBS.rotateZ(m, LIBS.degToRad(0));
                    LIBS.rotateY(m, LIBS.degToRad(30));
                    LIBS.rotateX(m, LIBS.degToRad(-100)); return m;
                })()
            },
            {
                parentName: 'leftHand', geom: Geometry.generateComplexCone(0.2, 0, 0.4, 4, 20, C.HEAD),
                trans: (() => {
                    let m = LIBS.get_I4();
                    // let m = createTransform(-1.225, 0.10, -0.425);
                    // LIBS.rotateZ(m, LIBS.degToRad(0));
                    // LIBS.rotateY(m, LIBS.degToRad(30));
                    LIBS.rotateX(m, LIBS.degToRad(-100));
                    return m;
                })()
            },
            {
                parentName: 'leftHand', geom: Geometry.generateComplexCone(0.2, 0, 0.4, 4, 20, C.HEAD),
                trans: (() => {
                    let m = LIBS.get_I4();
                    // let m = createTransform(-1.225, 0.16, -0.040);
                    // LIBS.rotateZ(m, LIBS.degToRad(180));
                    // LIBS.rotateY(m, LIBS.degToRad(30));
                    // LIBS.rotateX(m, LIBS.degToRad(-100));
                    return m;
                })()
            },

            {
                parentName: 'rightHand', geom: Geometry.generateHalfEllipsoid(0.15, 1.5, 0.6, 15, 15, C.HEAD, 'lower'),
                trans: (() => {
                    let m = LIBS.get_I4();
                    // let m = createTransform(1.225, 0.20, 0.2);
                    // LIBS.rotateZ(m, LIBS.degToRad(30));
                    // LIBS.rotateX(m, LIBS.degToRad(-10));
                    return m;
                })()
            },
            {
                parentName: 'rightHand', geom: Geometry.generateComplexCone(0.2, 0, 0.4, 4, 20, C.HEAD), trans: (() => {
                    let m = createTransform(0, 0, 0.55);
                    LIBS.rotateZ(m, LIBS.degToRad(180));
                    LIBS.rotateY(m, LIBS.degToRad(60));
                    LIBS.rotateX(m, LIBS.degToRad(-100));
                    return m;
                })()
            },
            {
                parentName: 'rightHand', geom: Geometry.generateComplexCone(0.2, 0, 0.4, 4, 20, C.HEAD), trans: (() => {
                    let m = createTransform(0, 0, 0.10);
                    LIBS.rotateZ(m, LIBS.degToRad(0));
                    LIBS.rotateY(m, LIBS.degToRad(60));
                    LIBS.rotateX(m, LIBS.degToRad(-100)); return m; })()
            },
            {
                parentName: 'rightHand', geom: Geometry.generateComplexCone(0.2, 0, 0.4, 4, 20, C.HEAD), trans: (() => {
                    let m = createTransform(0, 0, -0.55);
                    LIBS.rotateZ(m, LIBS.degToRad(0));
                    LIBS.rotateY(m, LIBS.degToRad(60));
                    LIBS.rotateX(m, LIBS.degToRad(-100)); return m; })()
            },
            {
                parentName: 'rightHand', geom: Geometry.generateComplexCone(0.2, 0, 0.4, 4, 20, C.HEAD), trans: (() => {
                    let m = createTransform(0, 0, 0.55);
                    LIBS.rotateZ(m, LIBS.degToRad(180));
                    LIBS.rotateY(m, LIBS.degToRad(60));
                    LIBS.rotateX(m, LIBS.degToRad(-100)); return m; })()
            },

            // Small fingers
            { parentName: 'leftHand', geom: Geometry.generateComplexCone(0.05, 0, 0.2, 50, 20, C.BEAK), trans: (() => {
                let m = createTransform(0.2, -0.9, 0);
                LIBS.rotateZ(m, LIBS.degToRad(-120));
                return m; })()},
            { parentName: 'leftHand', geom: Geometry.generateComplexCone(0.05, 0, 0.2, 50, 20, C.BEAK), trans: (() => {
                let m = createTransform(0.25, -0.8, 0.10);
                LIBS.rotateZ(m, LIBS.degToRad(-120)); return m; })()},
            { parentName: 'leftHand', geom: Geometry.generateComplexCone(0.05, 0, 0.2, 50, 20, C.BEAK), trans: (() => {
                let m = createTransform(0.22, -0.8, -0.10);
                LIBS.rotateZ(m, LIBS.degToRad(-120)); return m; })()},
            { parentName: 'rightHand', geom: Geometry.generateComplexCone(0.05, 0, 0.2, 50, 20, C.BEAK), trans: (() => {
                let m = createTransform(-0.2, -0.9, 0);
                LIBS.rotateZ(m, LIBS.degToRad(120)); return m; })()},
            { parentName: 'rightHand', geom: Geometry.generateComplexCone(0.05, 0, 0.2, 50, 20, C.BEAK), trans: (() => {
                let m = createTransform(-0.25, -0.8, 0.10);
                LIBS.rotateZ(m, LIBS.degToRad(120)); return m; })()},
            { parentName: 'rightHand', geom: Geometry.generateComplexCone(0.05, 0, 0.2, 50, 20, C.BEAK), trans: (() => {
                let m = createTransform(-0.22, -0.8, -0.10);
                LIBS.rotateZ(m, LIBS.degToRad(120)); return m; })()},


            // Legs
            { geom: Geometry.generateHalfHyperboloid(0.25, 0.5, 0.25, 20, 20, C.EMPO_LOWER_BODY, 0, 1.3), trans: (() => {
                    let m = createTransform(-0.55, -1.9, 0.3);
                    LIBS.rotateZ(m, LIBS.degToRad(-10));
                    LIBS.rotateX(m, LIBS.degToRad(-10));
                    return m;
                })()},

            { geom: Geometry.generateHalfHyperboloid(0.25, 0.5, 0.25, 20, 20, C.EMPO_LOWER_BODY, 0, 1.3), trans: (() => {
                    let m = createTransform(0.55, -1.9, 0.3);
                    LIBS.rotateZ(m, LIBS.degToRad(10));
                    LIBS.rotateX(m, LIBS.degToRad(-10));
                    return m;
                })()},

            // Feet
            { geom: Geometry.generateSphere(0.3, 0.12, 0.35, 10, 10, C.FEET), trans: createTransform(-0.55, -1.9, 0.24), animationType: 'none'},
            { geom: Geometry.generateSphere(0.3, 0.1, 0.5, 10, 10, C.FEET), trans: createTransform(-0.55, -2, 0.4), animationType: 'none'},

            { geom: Geometry.generateSphere(0.3, 0.12, 0.35, 10, 10, C.FEET), trans: createTransform(0.55, -1.9, 0.24), animationType: 'none'},
            { geom: Geometry.generateSphere(0.3, 0.1, 0.5, 10, 10, C.FEET), trans: createTransform(0.55, -2, 0.4), animationType: 'none'},

        ];

        // --- HIERARCHY BUILDING (using two-pass approach like Prinplup) ---
        const nodeMap = {}; // Stores all nodes by name

        // Pass 1: Create all ModelNode objects defined in partDefinitions
        partDefinitions.forEach(def => {
            // Give a default name if missing (useful for parts we don't animate/reference)
            if (!def.name) {
                def.name = `part_${Math.random().toString(16).substring(2)}`;
            }
            const node = new ModelNode(gl, def.geom, def.texture);
            if (def.trans) {
                node.setBaseTransform(def.trans);
            }
            nodeMap[def.name] = node;

            // Store references for animation
            if (def.name === 'leftHand') {
                this.animatedNodes.leftHand = node;
                this.baseTransforms.leftHand = def.trans; // Store base pose
            } else if (def.name === 'rightHand') {
                this.animatedNodes.rightHand = node;
                this.baseTransforms.rightHand = def.trans; // Store base pose
            }
            // Add other nodes to animatedNodes if needed
        });

        // Pass 2: Build the hierarchy based on parentName
        partDefinitions.forEach(def => {
            const node = nodeMap[def.name];
            const parentName = def.parentName || 'breathingNode'; // Default parent is breathingNode

            if (parentName === 'breathingNode') {
                breathingNode.addChild(node);
            } else {
                const parentNode = nodeMap[parentName];
                if (parentNode) {
                    parentNode.addChild(node);
                } else {
                    console.warn(`Parent node "${parentName}" not found for child "${def.name}". Attaching to breathingNode.`);
                    breathingNode.addChild(node); // Fallback parent
                }
            }
        });
        // --- END HIERARCHY BUILDING ---
    }


    draw(shader, parentMatrix) {
        // Combine the parent's matrix (e.g., mouse rotation)
        // with this model's static offset matrix
        const finalParentMatrix = LIBS.multiply(parentMatrix, this.modelMatrix);

        // Update all matrices in the tree
        this.rootNode.updateWorldMatrix(finalParentMatrix);

        // Start the recursive draw
        this.rootNode.draw(shader);
    }

    /**
     * Applies animations to the model's nodes.
     * @param {object} animValues - An object containing animation values (e.g., { body, flapAngle })
     */
    updateAnimation(animValues) {
        // Get the Y-axis animation value, defaulting to 0.0 if not provided
        const bodyY = animValues.body || 0.0;

        // Create a new translation matrix for the main body/breathing node
        const T_breath = LIBS.get_I4();

        // Apply the animation value PLUS the 2.5 static offset
        LIBS.translateY(T_breath, bodyY + 2.5);

        // Set this new matrix as the local transform for the *entire* model
        if (this.animatedNodes.breathingNode) {
            this.animatedNodes.breathingNode.setLocalTransform(T_breath);
        }

        // --- ADDED: Hand Flapping Animation ---
        const flapAngle = animValues.flapAngle || 0.0;
        // Adjust pivot points for Empoleon's hands (Y radius=1.5, X radius=0.2)
        const pivotY = 1.3; // Pivot near the top of the hand sphere
        const pivotX = 0.15; // Pivot near the edge of the hand sphere

        // --- Left Hand ---
        if (this.animatedNodes.leftHand && this.baseTransforms.leftHand) {
            // Create pivot matrices
            let T_up_left = LIBS.get_I4();
            LIBS.translateY(T_up_left, -pivotY);
            LIBS.translateX(T_up_left, pivotX);
            let T_down_left = LIBS.get_I4();
            LIBS.translateY(T_down_left, pivotY);
            LIBS.translateX(T_down_left, -pivotX);

            const R_left = LIBS.get_I4();
            LIBS.rotateY(R_left, flapAngle); // Rotate around Z axis

            // Combine for animation matrix: M_anim = T_up * R_z * T_down
            let leftAnim = LIBS.multiply(R_left, T_down_left);
            leftAnim = LIBS.multiply(T_up_left, leftAnim);

            // Combine with base: M_final = M_base * M_anim
            // IMPORTANT: Apply animation relative to the breathing node's transform
            // We apply the animation FIRST, then the base pose.
            const leftFinal = LIBS.multiply(this.baseTransforms.leftHand, leftAnim);
            this.animatedNodes.leftHand.setLocalTransform(leftFinal);
        }

        // --- Right Hand ---
        if (this.animatedNodes.rightHand && this.baseTransforms.rightHand) {
            // Create pivot matrices (mirrored X)
            let T_up_right = LIBS.get_I4();
            LIBS.translateY(T_up_right, -pivotY);
            LIBS.translateX(T_up_right, -pivotX); // Mirrored X
            let T_down_right = LIBS.get_I4();
            LIBS.translateY(T_down_right, pivotY);
            LIBS.translateX(T_down_right, pivotX); // Mirrored X

            const R_right = LIBS.get_I4();
            LIBS.rotateY(R_right, flapAngle); // Opposite direction

            // Combine for animation matrix: M_anim = T_up * R_z * T_down
            let rightAnim = LIBS.multiply(R_right, T_down_right);
            rightAnim = LIBS.multiply(T_up_right, rightAnim);

            // Combine with base: M_final = M_base * M_anim
            const rightFinal = LIBS.multiply(this.baseTransforms.rightHand, rightAnim);
            this.animatedNodes.rightHand.setLocalTransform(rightFinal);
        }
        // --- END ADDED ---
    }
}