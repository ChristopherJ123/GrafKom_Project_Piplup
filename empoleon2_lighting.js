// --- UTILITY FOR GEOMETRY ---
// We place the geometry generation logic in its own object to keep things organized.
const Geometry = {
    generateSphere: function (a, b, c, stack, step, color, textureOffset = 0) {
        const vertices = [];
        const faces = [];
        for (let i = 0; i <= stack; i++) {
            for (let j = 0; j <= step; j++) {
                const u = i / stack;
                const v = j / step;
                const theta = u * Math.PI;
                const phi = v * 2 * Math.PI;

                const x = a * Math.sin(theta) * Math.cos(phi);
                const y = b * Math.cos(theta);
                const z = c * Math.sin(theta) * Math.sin(phi);

                // Normal for a sphere is just its position vector (normalized)
                const normal = [x, y, z];

                vertices.push(x, y, z, color[0], color[1], color[2], v + textureOffset, u - 1.0, normal[0], normal[1], normal[2]);
            }
        }
        for (let i = 0; i < stack; i++) {
            for (let j = 0; j < step; j++) {
                const p1 = i * (step + 1) + j;
                const p2 = p1 + 1;
                const p3 = p1 + (step + 1);
                const p4 = p3 + 1;
                faces.push(p1, p2, p4, p1, p4, p3);
            }
        }
        return { vertices, faces };
    },

    generateTriBeak: function (baseWidth, height, depth, color) {
        const vertices = [];
        const faces = [];
        const halfBase = baseWidth / 2;

        // --- Define the 4 corner points ---
        const v0 = [-halfBase, -height / 2, 0];     // left bottom back
        const v1 = [halfBase, -height / 2, 0];      // right bottom back
        const v2 = [0, -height / 2, -depth];        // center bottom front
        const vTop = [0, height / 2, -depth / 2];   // top point

        // --- Calculate Normals for each Face ---

        // Face (0, 2, 3) [left]
        let A = [v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]]; // v2 - v0
        let B = [vTop[0] - v0[0], vTop[1] - v0[1], vTop[2] - v0[2]]; // vTop - v0
        let n_left = [
            A[1] * B[2] - A[2] * B[1],
            A[2] * B[0] - A[0] * B[2],
            A[0] * B[1] - A[1] * B[0]
        ];
        // Normalize
        let mag = Math.sqrt(n_left[0]**2 + n_left[1]**2 + n_left[2]**2);
        if (mag > 0) { n_left[0] /= mag; n_left[1] /= mag; n_left[2] /= mag; }

        // Face (1, 3, 2) [right]
        A = [vTop[0] - v1[0], vTop[1] - v1[1], vTop[2] - v1[2]]; // vTop - v1
        B = [v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2]]; // v2 - v1
        let n_right = [
            A[1] * B[2] - A[2] * B[1],
            A[2] * B[0] - A[0] * B[2],
            A[0] * B[1] - A[1] * B[0]
        ];
        // Normalize
        mag = Math.sqrt(n_right[0]**2 + n_right[1]**2 + n_right[2]**2);
        if (mag > 0) { n_right[0] /= mag; n_right[1] /= mag; n_right[2] /= mag; }

        // Face (0, 1, 2) [bottom]
        const n_bottom = [0, -1, 0]; // Points straight down

        // --- Push Vertices (3 per face) ---
        // Stride: x,y,z, r,g,b, u,v, nx,ny,nz (11 components)

        // Left Face (v0, v2, vTop)
        vertices.push(v0[0], v0[1], v0[2], color[0], color[1], color[2], 0, 0, n_left[0], n_left[1], n_left[2]);
        vertices.push(v2[0], v2[1], v2[2], color[0], color[1], color[2], 0, 0, n_left[0], n_left[1], n_left[2]);
        vertices.push(vTop[0], vTop[1], vTop[2], color[0], color[1], color[2], 0, 0, n_left[0], n_left[1], n_left[2]);

        // Right Face (v1, vTop, v2)
        vertices.push(v1[0], v1[1], v1[2], color[0], color[1], color[2], 0, 0, n_right[0], n_right[1], n_right[2]);
        vertices.push(vTop[0], vTop[1], vTop[2], color[0], color[1], color[2], 0, 0, n_right[0], n_right[1], n_right[2]);
        vertices.push(v2[0], v2[1], v2[2], color[0], color[1], color[2], 0, 0, n_right[0], n_right[1], n_right[2]);

        // Bottom Face (v0, v1, v2)
        vertices.push(v0[0], v0[1], v0[2], color[0], color[1], color[2], 0, 0, n_bottom[0], n_bottom[1], n_bottom[2]);
        vertices.push(v1[0], v1[1], v1[2], color[0], color[1], color[2], 0, 0, n_bottom[0], n_bottom[1], n_bottom[2]);
        vertices.push(v2[0], v2[1], v2[2], color[0], color[1], color[2], 0, 0, n_bottom[0], n_bottom[1], n_bottom[2]);

        // --- Push Faces ---
        // Each face just points to its 3 dedicated vertices
        faces.push(0, 1, 2); // Left
        faces.push(3, 4, 5); // Right
        faces.push(6, 7, 8); // Bottom

        return { vertices, faces };
    },

    generateTriangularPrism: function (baseWidth, height, depth, color) {
        const vertices = [];
        const faces = [];

        const halfBase = baseWidth / 2;
        const halfDepth = depth / 2;

        // --- Define the 6 corner points ---
        const f0 = [-halfBase, -height / 2, halfDepth]; // left bottom front
        const f1 = [halfBase, -height / 2, halfDepth];  // right bottom front
        const f2 = [0, height / 2, halfDepth];          // top front
        const b0 = [-halfBase, -height / 2, -halfDepth]; // left bottom back
        const b1 = [halfBase, -height / 2, -halfDepth];  // right bottom back
        const b2 = [0, height / 2, -halfDepth];          // top back

        // --- Define Normals ---
        const n_front = [0, 0, 1];
        const n_back = [0, 0, -1];
        const n_bottom = [0, -1, 0];

        // Normal for left face (f0, f2, b2, b0)
        let A = [f2[0] - f0[0], f2[1] - f0[1], f2[2] - f0[2]]; // f2 - f0
        let B = [b0[0] - f0[0], b0[1] - f0[1], b0[2] - f0[2]]; // b0 - f0
        let n_left = [
            A[1] * B[2] - A[2] * B[1],
            A[2] * B[0] - A[0] * B[2],
            A[0] * B[1] - A[1] * B[0]
        ];
        let mag = Math.sqrt(n_left[0]**2 + n_left[1]**2 + n_left[2]**2);
        if (mag > 0) { n_left[0] /= mag; n_left[1] /= mag; n_left[2] /= mag; }

        // Normal for right face (f1, b1, b2, f2)
        A = [b1[0] - f1[0], b1[1] - f1[1], b1[2] - f1[2]]; // b1 - f1
        B = [f2[0] - f1[0], f2[1] - f1[1], f2[2] - f1[2]]; // f2 - f1
        let n_right = [
            A[1] * B[2] - A[2] * B[1],
            A[2] * B[0] - A[0] * B[2],
            A[0] * B[1] - A[1] * B[0]
        ];
        mag = Math.sqrt(n_right[0]**2 + n_right[1]**2 + n_right[2]**2);
        if (mag > 0) { n_right[0] /= mag; n_right[1] /= mag; n_right[2] /= mag; }

        // Helper function to push 11-component vertex
        const pushVert = (v, n) => {
            vertices.push(v[0], v[1], v[2], color[0], color[1], color[2], 0, 0, n[0], n[1], n[2]);
        };

        // --- Create Vertices and Faces ---
        let v_idx = 0;

        // Front Face (f0, f1, f2)
        pushVert(f0, n_front);
        pushVert(f1, n_front);
        pushVert(f2, n_front);
        faces.push(v_idx, v_idx+1, v_idx+2);
        v_idx += 3;

        // Back Face (b0, b2, b1) - reverse winding
        pushVert(b0, n_back);
        pushVert(b1, n_back);
        pushVert(b2, n_back);
        faces.push(v_idx, v_idx+2, v_idx+1); // reversed
        v_idx += 3;

        // Bottom Face (f0, b0, b1, f1)
        pushVert(f0, n_bottom);
        pushVert(b0, n_bottom);
        pushVert(b1, n_bottom);
        pushVert(f1, n_bottom);
        faces.push(v_idx, v_idx+1, v_idx+2);
        faces.push(v_idx, v_idx+2, v_idx+3);
        v_idx += 4;

        // Left Face (f0, f2, b2, b0)
        pushVert(f0, n_left);
        pushVert(f2, n_left);
        pushVert(b2, n_left);
        pushVert(b0, n_left);
        faces.push(v_idx, v_idx+1, v_idx+2);
        faces.push(v_idx, v_idx+2, v_idx+3);
        v_idx += 4;

        // Right Face (f1, b1, b2, f2)
        pushVert(f1, n_right);
        pushVert(b1, n_right);
        pushVert(b2, n_right);
        pushVert(f2, n_right);
        faces.push(v_idx, v_idx+1, v_idx+2);
        faces.push(v_idx, v_idx+2, v_idx+3);
        // v_idx += 4; // no need

        return { vertices, faces };
    },



    generateHalfEllipsoid: function (a, b, c, stack, step, color, hemisphere = 'upper') {
        const vertices = [];
        const faces = [];

        for (let i = 0; i <= stack; i++) {
            for (let j = 0; j <= step; j++) {
                const u = i / stack;
                const v = j / step;

                let theta;
                if (hemisphere === 'upper') {
                    theta = u * (Math.PI / 2); // 0 to π/2
                } else if (hemisphere === 'lower') {
                    theta = (Math.PI / 2) + u * (Math.PI / 2); // π/2 to π
                } else {
                    theta = u * Math.PI; // full ellipsoid
                }

                const phi = v * 2 * Math.PI;

                const x = a * Math.sin(theta) * Math.cos(phi);
                const y = b * Math.cos(theta);
                const z = c * Math.sin(theta) * Math.sin(phi);

                // --- NEW: Calculate Normal ---
                let nx = x / (a * a);
                let ny = y / (b * b);
                let nz = z / (c * c);
                // Normalize the normal vector
                const mag = Math.sqrt(nx*nx + ny*ny + nz*nz);
                if (mag > 0) {
                    nx /= mag;
                    ny /= mag;
                    nz /= mag;
                }
                // --- End Normal Calculation ---

                // Push all 11 components
                vertices.push(x, y, z, color[0], color[1], color[2], v, u, nx, ny, nz);
            }
        }

        // Generate faces (unchanged)
        for (let i = 0; i < stack; i++) {
            for (let j = 0; j < step; j++) {
                const p1 = i * (step + 1) + j;
                const p2 = p1 + 1;
                const p3 = p1 + (step + 1);
                const p4 = p3 + 1;
                faces.push(p1, p2, p4, p1, p4, p3);
            }
        }

        return { vertices, faces };
    },

    generateComplexCone: function (radiusBottom, radiusTop, height, radialSegments, heightSegments, color, thetaStart = 0, thetaLength = Math.PI * 2) {
        const vertices = [];
        const faces = [];

        // --- NEW: Calculate the Y-component of the normal vector ---
        // This is based on the slope (change in radius vs. change in height)
        const deltaRadius = radiusBottom - radiusTop;
        const ny_component = deltaRadius;

        for (let i = 0; i <= heightSegments; i++) {
            for (let j = 0; j <= radialSegments; j++) {
                const u = i / heightSegments; // 0 to 1 (bottom to top)
                const v = j / radialSegments; // 0 to 1 (around circle)

                // Interpolate radius from bottom to top
                const radius = radiusBottom + (radiusTop - radiusBottom) * u;

                // Angle around the cone
                const theta = thetaStart + v * thetaLength;
                const cosTheta = Math.cos(theta);
                const sinTheta = Math.sin(theta);

                // Position
                const x = radius * cosTheta;
                const y = u * height - height / 2; // Center at origin
                const z = radius * sinTheta;

                // --- NEW: Calculate Normal ---
                // The normal's XZ components point radially from the Y-axis
                // The Y component is constant, based on the slope.
                let nx = height * cosTheta;
                let ny = ny_component;
                let nz = height * sinTheta;

                // Normalize the normal vector
                const mag = Math.sqrt(nx*nx + ny*ny + nz*nz);
                if (mag > 0) {
                    nx /= mag;
                    ny /= mag;
                    nz /= mag;
                }

                // Handle edge case of a flat disk
                if (height === 0) {
                    nx = 0;
                    ny = (radiusBottom > radiusTop) ? 1.0 : -1.0;
                    nz = 0;
                }
                // --- End Normal Calculation ---

                // Push all 11 components
                vertices.push(x, y, z, color[0], color[1], color[2], v, u, nx, ny, nz);
            }
        }

        // Generate faces (unchanged)
        for (let i = 0; i < heightSegments; i++) {
            for (let j = 0; j < radialSegments; j++) {
                const p1 = i * (radialSegments + 1) + j;
                const p2 = p1 + 1;
                const p3 = p1 + (radialSegments + 1);
                const p4 = p3 + 1;
                faces.push(p1, p2, p4, p1, p4, p3);
            }
        }

        return { vertices, faces };
    },

    generateHalfHyperboloid: function (a, b, c, stack, step, color, minHeight = 0, maxHeight = 1) {
        const vertices = [];
        const faces = [];

        for (let i = 0; i <= stack; i++) {
            for (let j = 0; j <= step; j++) {
                const u = (j / step) * 2 * Math.PI; // Full rotation
                const v = minHeight + (i / stack) * (maxHeight - minHeight); // Height

                const coshV = Math.cosh(v);
                const sinhV = Math.sinh(v);

                const x = a * coshV * Math.cos(u);
                const y = b * sinhV;
                const z = c * coshV * Math.sin(u);

                // --- NEW: Calculate Normal ---
                // Gradient is (2x/a², -2y/b², 2z/c²), which points inward.
                // We flip it to get the outward-facing normal.
                let nx = -x / (a * a);
                let ny = y / (b * b);
                let nz = -z / (c * c);
                // Normalize the normal vector
                const mag = Math.sqrt(nx*nx + ny*ny + nz*nz);
                if (mag > 0) {
                    nx /= mag;
                    ny /= mag;
                    nz /= mag;
                }
                // --- End Normal Calculation ---

                // Push all 11 components
                vertices.push(x, y, z, color[0], color[1], color[2], j / step, i / stack, nx, ny, nz);
            }
        }

        // Generate faces (unchanged)
        for (let i = 0; i < stack; i++) {
            for (let j = 0; j < step; j++) {
                const p1 = i * (step + 1) + j;
                const p2 = p1 + 1;
                const p3 = p1 + (step + 1);
                const p4 = p3 + 1;
                faces.push(p1, p2, p4, p1, p4, p3);
            }
        }

        return { vertices, faces };
    },


    generateBeak: function(width, thickness, length, segments, color) {
        const vertices = [];
        const faces = [];

        // Tip of the beak
        // Normal points straight out the tip
        vertices.push(0, 0, 0, color[0], color[1], color[2], 0, 0, 0, 0, 1);

        // Build the beak with circular cross-sections
        for (let i = 1; i <= segments; i++) {
            const t = i / segments; // Parameter from 0 to 1

            // Use sqrt(t) to make the beak fatter at the base and blunter at the tip
            const radiusScale = Math.sqrt(t);
            const currentZ = -length * t; // Move along the negative Z-axis

            for (let j = 0; j < segments; j++) {
                const theta = (j / segments) * 2 * Math.PI;
                const x = width * radiusScale * Math.cos(theta);
                const y = thickness * radiusScale * Math.sin(theta);

                // Normal for the sides
                let nx = x;
                let ny = y;
                const mag = Math.sqrt(nx * nx + ny * ny);
                if (mag > 0) {
                    nx /= mag;
                    ny /= mag;
                }

                vertices.push(x, y, currentZ, color[0], color[1], color[2], t, j / segments, nx, ny, 0);
            }
        }

        // Create faces for the tip
        for (let j = 1; j <= segments; j++) {
            faces.push(0, j, (j % segments) + 1);
        }

        // Create faces for the sides
        for (let i = 0; i < segments - 1; i++) {
            const ring1_start = 1 + i * segments;
            const ring2_start = 1 + (i + 1) * segments;
            for (let j = 0; j < segments; j++) {
                const p1 = ring1_start + j;
                const p2 = ring1_start + ((j + 1) % segments);
                const p3 = ring2_start + j;
                const p4 = ring2_start + ((j + 1) % segments);
                faces.push(p1, p3, p2, p2, p3, p4);
            }
        }

        return { vertices, faces };
    },

    generateTubeFromSpline: function(controlPoints, segments, radius, radialSegments, color) {
        var vertices = [];
        var faces = [];
        var splinePoints = [];
        var tangents = [];

        var points = [];
        points.push(controlPoints[0]);
        controlPoints.forEach(p => points.push(p));
        points.push(controlPoints[controlPoints.length - 1]);

        for (var i = 1; i < points.length - 2; i++) {
            var p0 = points[i - 1];
            var p1 = points[i];
            var p2 = points[i + 1];
            var p3 = points[i + 2];

            for (var j = 0; j <= segments; j++) {
                var t = j / segments;
                var t2 = t * t;
                var t3 = t2 * t;

                var x = 0.5 * ((2 * p1[0]) + (-p0[0] + p2[0]) * t + (2 * p0[0] - 5 * p1[0] + 4 * p2[0] - p3[0]) * t2 + (-p0[0] + 3 * p1[0] - 3 * p2[0] + p3[0]) * t3);
                var y = 0.5 * ((2 * p1[1]) + (-p0[1] + p2[1]) * t + (2 * p0[1] - 5 * p1[1] + 4 * p2[1] - p3[1]) * t2 + (-p0[1] + 3 * p1[1] - 3 * p2[1] + p3[1]) * t3);
                var z = 0.5 * ((2 * p1[2]) + (-p0[2] + p2[2]) * t + (2 * p0[2] - 5 * p1[2] + 4 * p2[2] - p3[2]) * t2 + (-p0[2] + 3 * p1[2] - 3 * p2[2] + p3[2]) * t3);
                splinePoints.push([x, y, z]);

                var tx = 0.5 * ((-p0[0] + p2[0]) + 2 * (2 * p0[0] - 5 * p1[0] + 4 * p2[0] - p3[0]) * t + 3 * (-p0[0] + 3 * p1[0] - 3 * p2[0] + p3[0]) * t2);
                var ty = 0.5 * ((-p0[1] + p2[1]) + 2 * (2 * p0[1] - 5 * p1[1] + 4 * p2[1] - p3[1]) * t + 3 * (-p0[1] + 3 * p1[1] - 3 * p2[1] + p3[1]) * t2);
                var tz = 0.5 * ((-p0[2] + p2[2]) + 2 * (2 * p0[2] - 5 * p1[2] + 4 * p2[2] - p3[2]) * t + 3 * (-p0[2] + 3 * p1[2] - 3 * p2[2] + p3[2]) * t2);

                var mag = Math.sqrt(tx * tx + ty * ty + tz * tz);
                tangents.push([tx / mag, ty / mag, tz / mag]);
            }
        }

        var up = [0, 1, 0];
        for (var i = 0; i < splinePoints.length; i++) {
            var point = splinePoints[i];
            var tangent = tangents[i];

            if (Math.abs(tangent[1]) > 0.999) {
                up = [1, 0, 0];
            } else {
                up = [0, 1, 0];
            }

            var normal = [tangent[1] * up[2] - tangent[2] * up[1], tangent[2] * up[0] - tangent[0] * up[2], tangent[0] * up[1] - tangent[1] * up[0]];
            var magN = Math.sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
            normal = [normal[0] / magN, normal[1] / magN, normal[2] / magN];

            var binormal = [tangent[1] * normal[2] - tangent[2] * normal[1], tangent[2] * normal[0] - tangent[0] * normal[2], tangent[0] * normal[1] - tangent[1] * normal[0]];

            for (var j = 0; j <= radialSegments; j++) {
                var theta = (j / radialSegments) * 2 * Math.PI;
                var x = point[0] + radius * (Math.cos(theta) * normal[0] + Math.sin(theta) * binormal[0]);
                var y = point[1] + radius * (Math.cos(theta) * normal[1] + Math.sin(theta) * binormal[1]);
                var z = point[2] + radius * (Math.cos(theta) * normal[2] + Math.sin(theta) * binormal[2]);

                // Normal is the vector from the spline point to the vertex
                var n = [x - point[0], y - point[1], z - point[2]];
                var mag = Math.sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
                if (mag > 0) { n[0] /= mag; n[1] /= mag; n[2] /= mag; }

                vertices.push(x, y, z, color[0], color[1], color[2], j / radialSegments, i / splinePoints.length, n[0], n[1], n[2]);
            }
        }

        for (var i = 0; i < splinePoints.length - 1; i++) {
            for (var j = 0; j < radialSegments; j++) {
                var p1 = i * (radialSegments + 1) + j;
                var p2 = p1 + 1;
                var p3 = (i + 1) * (radialSegments + 1) + j;
                var p4 = p3 + 1;
                faces.push(p1, p2, p4);
                faces.push(p1, p4, p3);
            }
        }
        return { vertices, faces };
    },

    generateTaperedShapeFromSpline: function(controlPoints, segments, startRadii, endRadii, radialSegments, color) {
        var vertices = [];
        var faces = [];
        var splinePoints = [];
        var tangents = [];

        // --- Spline Calculation (same as your existing function) ---
        var points = [];
        points.push(controlPoints[0]);
        controlPoints.forEach(p => points.push(p));
        points.push(controlPoints[controlPoints.length - 1]);

        for (var i = 1; i < points.length - 2; i++) {
            var p0 = points[i - 1];
            var p1 = points[i];
            var p2 = points[i + 1];
            var p3 = points[i + 2];

            for (var j = 0; j <= segments; j++) {
                var t = j / segments;
                var t2 = t * t;
                var t3 = t2 * t;

                var x = 0.5 * ((2 * p1[0]) + (-p0[0] + p2[0]) * t + (2 * p0[0] - 5 * p1[0] + 4 * p2[0] - p3[0]) * t2 + (-p0[0] + 3 * p1[0] - 3 * p2[0] + p3[0]) * t3);
                var y = 0.5 * ((2 * p1[1]) + (-p0[1] + p2[1]) * t + (2 * p0[1] - 5 * p1[1] + 4 * p2[1] - p3[1]) * t2 + (-p0[1] + 3 * p1[1] - 3 * p2[1] + p3[1]) * t3);
                var z = 0.5 * ((2 * p1[2]) + (-p0[2] + p2[2]) * t + (2 * p0[2] - 5 * p1[2] + 4 * p2[2] - p3[2]) * t2 + (-p0[2] + 3 * p1[2] - 3 * p2[2] + p3[2]) * t3);
                splinePoints.push([x, y, z]);

                var tx = 0.5 * ((-p0[0] + p2[0]) + 2 * (2 * p0[0] - 5 * p1[0] + 4 * p2[0] - p3[0]) * t + 3 * (-p0[0] + 3 * p1[0] - 3 * p2[0] + p3[0]) * t2);
                var ty = 0.5 * ((-p0[1] + p2[1]) + 2 * (2 * p0[1] - 5 * p1[1] + 4 * p2[1] - p3[1]) * t + 3 * (-p0[1] + 3 * p1[1] - 3 * p2[1] + p3[1]) * t2);
                var tz = 0.5 * ((-p0[2] + p2[2]) + 2 * (2 * p0[2] - 5 * p1[2] + 4 * p2[2] - p3[2]) * t + 3 * (-p0[2] + 3 * p1[2] - 3 * p2[2] + p3[2]) * t2);

                var mag = Math.sqrt(tx * tx + ty * ty + tz * tz);
                tangents.push([tx / mag, ty / mag, tz / mag]);
            }
        }

        // --- Mesh Generation (modified for tapering) ---
        var up = [0, 1, 0];
        for (var i = 0; i < splinePoints.length; i++) {
            var t = i / (splinePoints.length - 1); // Interpolation factor (0 to 1)

            // Interpolate radii for tapering effect
            var currentRadiusX = startRadii[0] * (1 - t) + endRadii[0] * t;
            var currentRadiusY = startRadii[1] * (1 - t) + endRadii[1] * t;

            var point = splinePoints[i];
            var tangent = tangents[i];

            if (Math.abs(tangent[1]) > 0.999) { up = [1, 0, 0]; }
            else { up = [0, 1, 0]; }

            var normal = [tangent[1] * up[2] - tangent[2] * up[1], tangent[2] * up[0] - tangent[0] * up[2], tangent[0] * up[1] - tangent[1] * up[0]];
            var magN = Math.sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
            normal = [normal[0] / magN, normal[1] / magN, normal[2] / magN];

            var binormal = [tangent[1] * normal[2] - tangent[2] * normal[1], tangent[2] * normal[0] - tangent[0] * normal[2], tangent[0] * normal[1] - tangent[1] * normal[0]];

            for (var j = 0; j <= radialSegments; j++) {
                var theta = (j / radialSegments) * 2 * Math.PI;
                // Use different radii for X and Y to create a flattened (elliptical) shape
                var x = point[0] + (currentRadiusX * Math.cos(theta) * normal[0] + currentRadiusY * Math.sin(theta) * binormal[0]);
                var y = point[1] + (currentRadiusX * Math.cos(theta) * normal[1] + currentRadiusY * Math.sin(theta) * binormal[1]);
                var z = point[2] + (currentRadiusX * Math.cos(theta) * normal[2] + currentRadiusY * Math.sin(theta) * binormal[2]);

                // Normal is vector from spline point to vertex
                var n = [x - point[0], y - point[1], z - point[2]];
                var mag = Math.sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
                if (mag > 0) { n[0] /= mag; n[1] /= mag; n[2] /= mag; }

                vertices.push(x, y, z, color[0], color[1], color[2], j/radialSegments, i/splinePoints.length, n[0], n[1], n[2]);
            }
        }

        for (var i = 0; i < splinePoints.length - 1; i++) {
            for (var j = 0; j < radialSegments; j++) {
                var p1 = i * (radialSegments + 1) + j;
                var p2 = p1 + 1;
                var p3 = (i + 1) * (radialSegments + 1) + j;
                var p4 = p3 + 1;
                faces.push(p1, p2, p4);
                faces.push(p1, p4, p3);
            }
        }
        return { vertices, faces };
    },

    generateLathe: function(points, segments, color) {
        const vertices = [];
        const faces = [];
        const phi_step = 2 * Math.PI / segments;

        // Generate vertices
        for (let i = 0; i <= segments; i++) {
            const phi = i * phi_step;
            const sin_phi = Math.sin(phi);
            const cos_phi = Math.cos(phi);

            for (let j = 0; j < points.length; j++) {
                const x = points[j][0] * cos_phi;
                const y = points[j][1];
                const z = points[j][0] * sin_phi;

                // Calculate normal based on profile tangent
                let p_prev = points[j > 0 ? j - 1 : j];
                let p_next = points[j < points.length - 1 ? j + 1 : j];

                let tangent_r, tangent_y;
                if (j == 0) {
                    tangent_r = points[1][0] - points[0][0];
                    tangent_y = points[1][1] - points[0][1];
                } else if (j == points.length - 1) {
                    tangent_r = points[j][0] - points[j-1][0];
                    tangent_y = points[j][1] - points[j-1][1];
                } else {
                    tangent_r = p_next[0] - p_prev[0];
                    tangent_y = p_next[1] - p_prev[1];
                }

                let profile_normal_r = -tangent_y;
                let profile_normal_y = tangent_r;

                let mag = Math.sqrt(profile_normal_r*profile_normal_r + profile_normal_y*profile_normal_y);
                if (mag > 0) {
                    profile_normal_r /= mag;
                    profile_normal_y /= mag;
                }

                let nx = profile_normal_r * cos_phi;
                let ny = profile_normal_y;
                let nz = profile_normal_r * sin_phi;
                // End normal calculation

                // UV coordinates
                const u = i / segments;
                const v = j / (points.length - 1);

                vertices.push(x, y, z, color[0], color[1], color[2], u, v, nx, ny, nz);
            }
        }

        // Generate faces
        for (let i = 0; i < segments; i++) {
            for (let j = 0; j < points.length - 1; j++) {
                const p1 = i * points.length + j;
                const p2 = p1 + points.length;
                const p3 = p1 + 1;
                const p4 = p2 + 1;

                faces.push(p1, p2, p4, p1, p4, p3);
            }
        }

        return { vertices, faces };
    },
};

// --- PIPLUP PART CLASS ---
// Represents a single drawable part of the Piplup model.
class ModelNode {
    constructor(gl, geometry = null, texture = null) {
        this.gl = gl;
        this.geometry = geometry; // The drawable geometry
        this.texture = texture;
        this.buffers = null;

        this.baseMatrix = LIBS.get_I4();  // The "default" pose (set once)
        this.localMatrix = LIBS.get_I4(); // The final transform (base * animation)
        this.worldMatrix = LIBS.get_I4(); // Final transformation in world space
        this.children = [];
        this.parent = null;

        if (this.geometry) {
            this.buffers = this.createBuffers();
        }
    }

    addChild(node) {
        node.parent = this;
        this.children.push(node);
    }

    createBuffers() {
        const vertexBuffer = this.gl.createBuffer();
        this.gl.bindBuffer(this.gl.ARRAY_BUFFER, vertexBuffer);
        this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(this.geometry.vertices), this.gl.STATIC_DRAW);

        const facesBuffer = this.gl.createBuffer();
        this.gl.bindBuffer(this.gl.ELEMENT_ARRAY_BUFFER, facesBuffer);
        this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(this.geometry.faces), this.gl.STATIC_DRAW);

        return { vertex: vertexBuffer, faces: facesBuffer, faces_length: this.geometry.faces.length };
    }

    setBaseTransform(matrix) {
        this.baseMatrix = matrix;
        this.localMatrix = matrix; // By default, local = base
    }

    setLocalTransform(matrix) {
        this.localMatrix = matrix;
    }

    updateWorldMatrix(parentWorldMatrix) {
        if (parentWorldMatrix) {
            this.worldMatrix = LIBS.multiply(this.localMatrix, parentWorldMatrix);
        } else {
            this.worldMatrix = this.localMatrix;
        }

        for (const child of this.children) {
            child.updateWorldMatrix(this.worldMatrix);
        }
    }

    draw(shader) {
        if (this.buffers) {
            const gl = this.gl;
            gl.uniformMatrix4fv(shader.locations.Mmatrix, false, this.worldMatrix);
            gl.bindBuffer(gl.ARRAY_BUFFER, this.buffers.vertex);

            // **IMPORTANT: Stride is now 11 floats (3pos, 3color, 2uv, 3normal)**
            const stride = 4 * (3 + 3 + 2 + 3);
            gl.vertexAttribPointer(shader.locations.position, 3, gl.FLOAT, false, stride, 0);
            gl.vertexAttribPointer(shader.locations.color, 3, gl.FLOAT, false, stride, 3 * 4);
            gl.vertexAttribPointer(shader.locations.texcoord, 2, gl.FLOAT, false, stride, 6 * 4);

            // **NEW: Enable the normal attribute**
            if (shader.locations.normal !== -1) {
                gl.vertexAttribPointer(shader.locations.normal, 3, gl.FLOAT, false, stride, 8 * 4);
            }

            if (this.texture) {
                gl.activeTexture(gl.TEXTURE0);
                gl.bindTexture(gl.TEXTURE_2D, this.texture);
                gl.uniform1i(shader.locations.sampler, 0);
                gl.uniform1i(shader.locations.u_useTexture, 1);
            } else {
                gl.uniform1i(shader.locations.u_useTexture, 0);
            }

            gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, this.buffers.faces);
            gl.drawElements(gl.TRIANGLES, this.buffers.faces_length, gl.UNSIGNED_SHORT, 0);
        }

        for (const child of this.children) {
            child.draw(shader);
        }
    }
}
// --- Empoleon CONTAINER CLASS ---
// Manages all the parts that make up the Empoleon.
class Empoleon {
    constructor(gl, renderer) {
        this.gl = gl;
        this.renderer = renderer;

        // NEW: Root of the scene graph
        this.rootNode = new ModelNode(gl);

        // This will be its static offset matrix (like Y+2.5)
        this.modelMatrix = LIBS.get_I4();

        // NEW: For future animations
        this.animatedNodes = {};
        this.baseTransforms = {};

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

        // First, define the 2D profile for the bean shape
        // const body_profile = [
        //     [0.0, 1.9, 0],
        //     // [0.1, 1.95, 0],
        //     [0.25, 1.9, 0],
        //     [0.4, 1.8, 0],
        //     [0.5, 1.6, 0],
        //     [0.52, 1.4, 0],
        //     [0.54, 1.2, 0],   // Top point (neck)
        //     [0.6, 0.9, 0],
        //     [0.8, 0.4, 0],
        //     [0.95, -0.1, 0],   // Widest part of the belly
        //     [0.9, -0.5, 0],
        //     [0.6, -0.9, 0],
        //     // [0.45, -1.1, 0],
        //     [0.0, -1.1, 0]   // Bottom point
        // ];

        const body_profile = [
            [0.0, 1.95, 0],
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
            [0.0, -1.3, 0],
            [0.0, -1.3, 0],   // Bottom point
        ];

        // Define parts and their local transformations
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
                    // Control points define the curve's path from base to tip
                    [
                        [0.0, -0.6, -0.8],  // Start point on the lower back
                        [0.0, -0.6, -1.2],  // Mid-point, curving down and back
                        [0.0, -0.6, -2.0]   // End point, the tip of the tail
                    ],
                    50,   // Segments for a smooth curve
                    [0.2, 0.9],  // Start Radii [thickness, width] - wide and flat at the base
                    [0.01, 0.01], // End Radii [thickness, width] - narrow and thin at the tip
                    20,   // Radial segments
                    C.HEAD // Using the head color for the tail
                ),
                trans: (() => {
                    const m = createTransform(0, 0.6, 0.1);
                    // LIBS.scale(m, 1.2);
                    return m
                })(), // No transformation needed, points are in world space
            },

            {
                geom: Geometry.generateTaperedShapeFromSpline(
                    // Control points define the curve's path from base to tip
                    [
                        [0.0, 0.2, 1.0],  // Start point on the lower back
                        [0.0, 0.2, 1.2],  // Mid-point, curving down and back
                        [0.0, 0.2, 1.4]   // End point, the tip of the tail
                    ],
                    50,   // Segments for a smooth curve
                    [0.15, 0.9],  // Start Radii [thickness, width] - wide and flat at the base
                    [0.01, 0.01], // End Radii [thickness, width] - narrow and thin at the tip
                    20,   // Radial segments
                    C.HEAD // Using the head color for the tail
                ),
                trans: (() => {
                    const m = createTransform(0, -0.35, 0.1);
                    LIBS.rotateX(m, -Math.PI / 8);
                    // LIBS.scale(m, 1.2);
                    return m
                })(), // No transformation needed, points are in world space
            },

            // BEAK
            // { geom: Geometry.generateCone(0.3, 0, 0.5, 3, 20, C.BEAK), trans: (() => {
            //         let m = createTransform(0, 1.6, 0.55);
            //             LIBS.rotateZ(m, LIBS.degToRad(0));
            //             LIBS.rotateY(m, LIBS.degToRad(30));
            //             LIBS.rotateX(m, LIBS.degToRad(75));
            //             return m;
            //     })()},

            // { geom: Geometry.generateCone(0.4, 0, 0.7, 3, 20, C.BEAK), trans: (() => {
            //         let m = createTransform(0, 1.75, 0.45);
            //             LIBS.rotateZ(m, LIBS.degToRad(0));
            //             LIBS.rotateY(m, LIBS.degToRad(-30));
            //             LIBS.rotateX(m, LIBS.degToRad(105));
            //             return m;
            //     })()},

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
            // Tengah (head/body decor)
            { geom: Geometry.generateTriangularPrism(0.1, 0.1, 1.15, C.BEAK), trans: (() => {
                    let m = createTransform(0, 2.37, 0.45);
                        LIBS.rotateZ(m, LIBS.degToRad(180));
                        LIBS.rotateY(m, LIBS.degToRad(0));
                        LIBS.rotateX(m, LIBS.degToRad(-100));
                        return m;
                })()},

            // Kiri (head/body decor)
            { geom: Geometry.generateTriangularPrism(0.1, 0.1, 0.7, C.BEAK), trans: (() => {
                    let m = createTransform(0.43, 2.4, 0.20);
                        LIBS.rotateZ(m, LIBS.degToRad(210));
                        LIBS.rotateY(m, LIBS.degToRad(0));
                        LIBS.rotateX(m, LIBS.degToRad(-100));
                        return m;
                })()},

            // Kanan (head/body decor)
            { geom: Geometry.generateTriangularPrism(0.1, 0.1, 0.7, C.BEAK), trans: (() => {
                    let m = createTransform(-0.43, 2.4, 0.20);
                        LIBS.rotateZ(m, LIBS.degToRad(-210));
                        LIBS.rotateY(m, LIBS.degToRad(0));
                        LIBS.rotateX(m, LIBS.degToRad(-100));
                        return m;
                })()},

            // Sambungan Beak-Crown Kiri (head/body decor)
            { geom: Geometry.generateTriangularPrism(0.1, 0.1, 0.5, C.BEAK), trans: (() => {
                    let m = createTransform(0.30, 1.9, 0.37);
                        LIBS.rotateZ(m, LIBS.degToRad(-15));
                        LIBS.rotateY(m, LIBS.degToRad(33));
                        LIBS.rotateX(m, LIBS.degToRad(-125));
                        return m;
                })()},

            // Sambungan Beak-Crown Kanan (head/body decor)
            { geom: Geometry.generateTriangularPrism(0.1, 0.1, 0.5, C.BEAK), trans: (() => {
                    let m = createTransform(-0.30, 1.9, 0.37);
                        LIBS.rotateZ(m, LIBS.degToRad(15));
                        LIBS.rotateY(m, LIBS.degToRad(-33));
                        LIBS.rotateX(m, LIBS.degToRad(-125));
                        return m;
                })()},

            // Crown Diamond Tengah (head/body decor)
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

            // Crown Diamond Kanan (head/body decor)
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

            // Crown Diamond Kiri (head/body decor)
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

            // KERAH Kiri (head/body decor)
            {
                geom: Geometry.generateTaperedShapeFromSpline(
                    // Control points define the curve's path from base to tip
                    [
                        [0.0, 0.2, 1.4],  // Start point on the lower back
                        [0.0, 0.2, 1.2],  // Mid-point, curving down and back
                        [0.0, 0.2, 1.0]   // End point, the tip of the tail
                    ],
                    50,   // Segments for a smooth curve
                    [0.02, 0.5],  // Start Radii [thickness, width] - wide and flat at the base
                    [0.01, 0.01], // End Radii [thickness, width] - narrow and thin at the tip
                    20,   // Radial segments
                    C.HEAD // Using the head color for the tail
                ),
                trans: (() => {
                    const m = createTransform(0.17, 1.72 , -0.85);
                    LIBS.rotateX(m, LIBS.degToRad(-15));
                    LIBS.rotateY(m, LIBS.degToRad(-15));
                    LIBS.rotateZ(m, LIBS.degToRad(40));

                    // LIBS.scale(m, 1.2);
                    return m
                })(), // No transformation needed, points are in world space
            },

            {
                geom: Geometry.generateTaperedShapeFromSpline(
                    // Control points define the curve's path from base to tip
                    [
                        [0.0, 0.2, 1.4],  // Start point on the lower back
                        [0.0, 0.2, 1.2],  // Mid-point, curving down and back
                        [0.0, 0.2, 1.0]   // End point, the tip of the tail
                    ],
                    50,   // Segments for a smooth curve
                    [0.02, 0.5],  // Start Radii [thickness, width] - wide and flat at the base
                    [0.01, 0.01], // End Radii [thickness, width] - narrow and thin at the tip
                    20,   // Radial segments
                    C.HEAD // Using the head color for the tail
                ),
                trans: (() => {
                    const m = createTransform(-0.17, 1.72 , -0.85);
                    LIBS.rotateX(m, LIBS.degToRad(-15));
                    LIBS.rotateY(m, LIBS.degToRad(15));
                    LIBS.rotateZ(m, LIBS.degToRad(-40));

                    // LIBS.scale(m, 1.2);
                    return m
                })(), // No transformation needed, points are in world space
            },

            // {
            //     geom: Geometry.generateTaperedShapeFromSpline(
            //         // Control points define the curve's path from base to tip
            //         [
            //             [0.0, 0.2, 1.0],  // Start point on the lower back
            //             [0.0, 0.2, 1.2],  // Mid-point, curving down and back
            //             [0.0, 0.2, 1.4]   // End point, the tip of the tail
            //         ],
            //         50,   // Segments for a smooth curve
            //         [0.25, 0.6],  // Start Radii [thickness, width] - wide and flat at the base
            //         [0.01, 0.01], // End Radii [thickness, width] - narrow and thin at the tip
            //         20,   // Radial segments
            //         C.HEAD // Using the head color for the tail
            //     ),
            //     trans: (() => {
            //         const m = createTransform(0, -0.65, 0.3);
            //         LIBS.rotateX(m, -0.5);
            //         // LIBS.scale(m, 1.2);
            //         return m
            //     })(), // No transformation needed, points are in world space
            // },

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

            // Hands
            { geom: Geometry.generateSphere(0.2, 1.5, 0.5, 15, 15, C.EMPO_BASE), trans: (() => {
                    let m = createTransform(-1.2, 0.25, 0.2);
                    LIBS.rotateZ(m, LIBS.degToRad(-30));
                    LIBS.rotateX(m, LIBS.degToRad(-10));
                    return m;
                })()},
            { geom: Geometry.generateSphere(0.2, 1.5, 0.5, 15, 15, C.EMPO_BASE), trans: (() => {
                    let m = createTransform(1.2, 0.25, 0.2);
                    LIBS.rotateZ(m, LIBS.degToRad(30));
                    LIBS.rotateX(m, LIBS.degToRad(-10));
                    return m;
                })()},

            // Outer Hands
            // Kiri
            { geom: Geometry.generateHalfEllipsoid(0.15, 1.5, 0.6, 15, 15, C.HEAD, 'lower'), trans: (() => {
                    let m = createTransform(-1.225, 0.20, 0.2);
                    LIBS.rotateZ(m, LIBS.degToRad(-30));
                    LIBS.rotateX(m, LIBS.degToRad(-10));
                    return m;
                })()},

            { geom: Geometry.generateComplexCone(0.2, 0, 0.4, 4, 20, C.HEAD), trans: (() => {
                    let m = createTransform(-1.225, 0.30, 0.85);
                        LIBS.rotateZ(m, LIBS.degToRad(180));
                        LIBS.rotateY(m, LIBS.degToRad(30));
                        LIBS.rotateX(m, LIBS.degToRad(-100));
                        return m;
                })()},

            { geom: Geometry.generateComplexCone(0.2, 0, 0.4, 4, 20, C.HEAD), trans: (() => {
                    let m = createTransform(-1.225, 0.225, 0.465);
                        LIBS.rotateZ(m, LIBS.degToRad(0));
                        LIBS.rotateY(m, LIBS.degToRad(30));
                        LIBS.rotateX(m, LIBS.degToRad(-100));
                        return m;
                })()},

            { geom: Geometry.generateComplexCone(0.2, 0, 0.4, 4, 20, C.HEAD), trans: (() => {
                    let m = createTransform(-1.225, 0.10, -0.425);
                        LIBS.rotateZ(m, LIBS.degToRad(0));
                        LIBS.rotateY(m, LIBS.degToRad(30));
                        LIBS.rotateX(m, LIBS.degToRad(-100));
                        return m;
                })()},

            { geom: Geometry.generateComplexCone(0.2, 0, 0.4, 4, 20, C.HEAD), trans: (() => {
                    let m = createTransform(-1.225, 0.16, -0.040);
                        LIBS.rotateZ(m, LIBS.degToRad(180));
                        LIBS.rotateY(m, LIBS.degToRad(30));
                        LIBS.rotateX(m, LIBS.degToRad(-100));
                        return m;
                })()},

            // Kanan
            { geom: Geometry.generateHalfEllipsoid(0.15, 1.5, 0.6, 15, 15, C.HEAD, 'lower'), trans: (() => {
                    let m = createTransform(1.225, 0.20, 0.2);
                    LIBS.rotateZ(m, LIBS.degToRad(30));
                    LIBS.rotateX(m, LIBS.degToRad(-10));
                    return m;
                })()},

            { geom: Geometry.generateComplexCone(0.2, 0, 0.4, 4, 20, C.HEAD), trans: (() => {
                    let m = createTransform(1.225, 0.30, 0.85);
                        LIBS.rotateZ(m, LIBS.degToRad(180));
                        LIBS.rotateY(m, LIBS.degToRad(60));
                        LIBS.rotateX(m, LIBS.degToRad(-100));
                        return m;
                })()},

            { geom: Geometry.generateComplexCone(0.2, 0, 0.4, 4, 20, C.HEAD), trans: (() => {
                    let m = createTransform(1.225, 0.225, 0.465);
                        LIBS.rotateZ(m, LIBS.degToRad(0));
                        LIBS.rotateY(m, LIBS.degToRad(60));
                        LIBS.rotateX(m, LIBS.degToRad(-100));
                        return m;
                })()},

            { geom: Geometry.generateComplexCone(0.2, 0, 0.4, 4, 20, C.HEAD), trans: (() => {
                    let m = createTransform(1.225, 0.10, -0.425);
                        LIBS.rotateZ(m, LIBS.degToRad(0));
                        LIBS.rotateY(m, LIBS.degToRad(60));
                        LIBS.rotateX(m, LIBS.degToRad(-100));
                        return m;
                })()},

            { geom: Geometry.generateComplexCone(0.2, 0, 0.4, 4, 20, C.HEAD), trans: (() => {
                    let m = createTransform(1.22, 0.16, -0.040);
                        LIBS.rotateZ(m, LIBS.degToRad(180));
                        LIBS.rotateY(m, LIBS.degToRad(60));
                        LIBS.rotateX(m, LIBS.degToRad(-100));
                        return m;
                })()},

            // Finger Kanan
            { geom: Geometry.generateComplexCone(0.05, 0, 0.2, 50, 20, C.BEAK), trans: (() => {
                    let m = createTransform(-1.7, -0.9, 0.40);
                        LIBS.rotateZ(m, LIBS.degToRad(-120));
                        LIBS.rotateY(m, LIBS.degToRad(0));
                        LIBS.rotateX(m, LIBS.degToRad(0));
                        return m;
                })()},

            // Finger Kanan
            { geom: Geometry.generateComplexCone(0.05, 0, 0.2, 50, 20, C.BEAK), trans: (() => {
                    let m = createTransform(-1.65, -0.8, 0.50);
                        LIBS.rotateZ(m, LIBS.degToRad(-120));
                        LIBS.rotateY(m, LIBS.degToRad(0));
                        LIBS.rotateX(m, LIBS.degToRad(0));
                        return m;
                })()},

            // Finger Kanan
            { geom: Geometry.generateComplexCone(0.05, 0, 0.2, 50, 20, C.BEAK), trans: (() => {
                    let m = createTransform(-1.62, -0.8, 0.30);
                        LIBS.rotateZ(m, LIBS.degToRad(-120));
                        LIBS.rotateY(m, LIBS.degToRad(0));
                        LIBS.rotateX(m, LIBS.degToRad(0));
                        return m;
                })()},

            // Finger Kiri
            { geom: Geometry.generateComplexCone(0.05, 0, 0.2, 50, 20, C.BEAK), trans: (() => {
                    let m = createTransform(1.7, -0.9, 0.40);
                        LIBS.rotateZ(m, LIBS.degToRad(120));
                        LIBS.rotateY(m, LIBS.degToRad(0));
                        LIBS.rotateX(m, LIBS.degToRad(0));
                        return m;
                })()},

            // Finger Kiri
            { geom: Geometry.generateComplexCone(0.05, 0, 0.2, 50, 20, C.BEAK), trans: (() => {
                    let m = createTransform(1.65, -0.8, 0.50);
                        LIBS.rotateZ(m, LIBS.degToRad(120));
                        LIBS.rotateY(m, LIBS.degToRad(0));
                        LIBS.rotateX(m, LIBS.degToRad(0));
                        return m;
                })()},

            // Finger Kiri
            { geom: Geometry.generateComplexCone(0.05, 0, 0.2, 50, 20, C.BEAK), trans: (() => {
                    let m = createTransform(1.62, -0.8, 0.30);
                        LIBS.rotateZ(m, LIBS.degToRad(120));
                        LIBS.rotateY(m, LIBS.degToRad(0));
                        LIBS.rotateX(m, LIBS.degToRad(0));
                        return m;
                })()},

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

    ]
        // --- NEW: Convert loop to use ModelNode ---
        partDefinitions.forEach(def => {
            // Create a ModelNode instead of a PiplupPart
            const node = new ModelNode(gl, def.geom, def.texture);

            // Set its base transform
            node.setBaseTransform(def.trans);

            // Add it as a child of the root
            this.rootNode.addChild(node);
        });
    }

    draw(shader, parentMatrix) {
        // Combine the parent's matrix (e.g., mouse rotation)
        // with this model's static offset matrix
        const finalParentMatrix = LIBS.multiply(parentMatrix, this.modelMatrix);

        // Update all matrices in the tree
        this.rootNode.updateWorldMatrix(finalParentMatrix);

        // Start the recursive draw
        this.rootNode.draw(shader);
    }}

// --- RENDERER CLASS ---
// Manages the overall WebGL scene, shaders, and render loop.
class Renderer {
    constructor(canvasId) {
        this.canvas = document.getElementById(canvasId);
        this.canvas.width = window.innerWidth;
        this.canvas.height = window.innerHeight;

        this.gl = this.canvas.getContext("webgl", { antialias: true });
        if (!this.gl) throw new Error("WebGL not supported");

        this.shader = this.createShaderProgram();
        this.empoleon = new Empoleon(this.gl, this);

        this.viewMatrix = LIBS.get_I4();
        LIBS.translateZ(this.viewMatrix, -9);
        // LIBS.rotateY(this.viewMatrix, Math.PI / 2);
        this.projMatrix = LIBS.get_projection(40, this.canvas.width / this.canvas.height, 1, 100);

        // New mouse rotation matrix
        this.mouseRotationMatrix = LIBS.get_I4();

        this.initInputHandlers();
        this.startRenderLoop();
    }

    createShaderProgram() {
        const gl = this.gl;
        const vsSource = `
            attribute vec3 position;
            attribute vec3 color;
            attribute vec2 texcoord;
            attribute vec3 normal; // <-- ADD THIS
            uniform mat4 Mmatrix, Vmatrix, Pmatrix;
            varying vec3 vColor;
            varying vec2 vTexcoord;
            void main(void) {
                gl_Position = Pmatrix * Vmatrix * Mmatrix * vec4(position, 1.);
                vColor = color;
                vTexcoord = texcoord;
            }`;
        const fsSource = `
            precision mediump float;
            varying vec3 vColor;
            varying vec2 vTexcoord;
            uniform sampler2D sampler;
            uniform int u_useTexture;

            void main(void) {
                if (u_useTexture == 1) {
                    gl_FragColor = texture2D(sampler, vTexcoord);
                } else {
                    gl_FragColor = vec4(vColor, 1.);
                }
            }`;

        const vs = this.compileShader(vsSource, gl.VERTEX_SHADER);
        const fs = this.compileShader(fsSource, gl.FRAGMENT_SHADER);

        const program = gl.createProgram();
        gl.attachShader(program, vs);
        gl.attachShader(program, fs);
        gl.linkProgram(program);

        gl.useProgram(program);

        const locations = {
            position: gl.getAttribLocation(program, "position"),
            color: gl.getAttribLocation(program, "color"),
            texcoord: gl.getAttribLocation(program, "texcoord"),
            normal: gl.getAttribLocation(program, "normal"),
            Pmatrix: gl.getUniformLocation(program, "Pmatrix"),
            Vmatrix: gl.getUniformLocation(program, "Vmatrix"),
            Mmatrix: gl.getUniformLocation(program, "Mmatrix"),
            sampler: gl.getUniformLocation(program, "sampler"),
            u_useTexture: gl.getUniformLocation(program, "u_useTexture")
        };

        gl.enableVertexAttribArray(locations.position);
        gl.enableVertexAttribArray(locations.color);
        gl.enableVertexAttribArray(locations.texcoord);
        gl.enableVertexAttribArray(locations.normal);

        return { program, locations };
    }

    compileShader(source, type) {
        const gl = this.gl;
        const shader = gl.createShader(type);
        gl.shaderSource(shader, source);
        gl.compileShader(shader);
        if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
            throw new Error("Shader compile error: " + gl.getShaderInfoLog(shader));
        }
        return shader;
    }

    loadTexture(url) {
        const gl = this.gl;
        const texture = gl.createTexture();
        gl.bindTexture(gl.TEXTURE_2D, texture);

        // Placeholder pixel
        gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, 1, 1, 0, gl.RGBA, gl.UNSIGNED_BYTE, new Uint8Array([0, 0, 255, 255]));

        const image = new Image();
        image.onload = function () {
            gl.bindTexture(gl.TEXTURE_2D, texture);
            gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, image);
            gl.generateMipmap(gl.TEXTURE_2D);
        };
        image.src = url;

        return texture;
    }

    initInputHandlers() {
        let drag = false;
        let x_prev, y_prev;
        let dX = 0, dY = 0;
        let THETA = 0, PHI = 0;
        const FRICTION = 0.15;

        this.canvas.onmousedown = (e) => { drag = true; x_prev = e.pageX; y_prev = e.pageY; };
        this.canvas.onmouseup = () => { drag = false; };
        this.canvas.onmouseout = () => { drag = false; };
        this.canvas.onmousemove = (e) => {
            if (!drag) return;
            dX = (e.pageX - x_prev) * 2 * Math.PI / this.canvas.width;
            dY = (e.pageY - y_prev) * 2 * Math.PI / this.canvas.height;
            THETA += dX;
            PHI += dY;
            x_prev = e.pageX;
            y_prev = e.pageY;
        };

        // This function is called every frame to update the Piplup's rotation
        this.updateRotation = () => {
            if (!drag) {
                dX *= (1 - FRICTION);
                dY *= (1 - FRICTION);
                THETA += dX;
                PHI += dY;
            }
            const rotationMatrix = LIBS.get_I4();
            LIBS.rotateY(rotationMatrix, THETA);
            LIBS.rotateX(rotationMatrix, PHI);
            this.mouseRotationMatrix = rotationMatrix;
        };
    }

    startRenderLoop() {
        const gl = this.gl;
        gl.enable(gl.DEPTH_TEST);
        gl.depthFunc(gl.LEQUAL);
        gl.clearColor(0.0, 0.0, 0.0, 0.0);
        gl.clearDepth(1.0);

        const render = () => {
            this.updateRotation();

            gl.viewport(0, 0, this.canvas.width, this.canvas.height);
            gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

            gl.uniformMatrix4fv(this.shader.locations.Pmatrix, false, this.projMatrix);
            gl.uniformMatrix4fv(this.shader.locations.Vmatrix, false, this.viewMatrix);

            this.empoleon.draw(this.shader, this.mouseRotationMatrix);

            requestAnimationFrame(render);
        };
        render();
    }
}

// --- START THE APPLICATION ---
window.addEventListener('load', () => {
    new Renderer('myCanvas');
});