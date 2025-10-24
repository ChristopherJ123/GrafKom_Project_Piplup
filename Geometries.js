const Geometries = {
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

    /**
     * Generates a blunted, beak-like shape.
     * @param {number} width - Max width (x-axis).
     * @param {number} thickness - Max thickness (y-axis).
     * @param {number} length - The length of the beak (z-axis).
     * @param {number} segments - The number of segments for resolution.
     * @param {Array<number>} color - The RGB color array.
     * @returns {{vertices: Array<number>, faces: Array<number>}}
     */
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

    /**
     * Generates one half (the top half, in positive Y) of a blunted, beak-like shape,
     * including vertex normals for lighting.
     * @param {number} width - Max width (x-axis).
     * @param {number} thickness - Max thickness (y-axis).
     * @param {number} length - The length of the beak (z-axis).
     * @param {number} segments - The number of segments for resolution.
     * @param {Array<number>} color - The RGB color array.
     * @returns {{vertices: Array<number>, faces: Array<number>}}
     */
    generateBeakHalf: function(width, thickness, length, segments, color) {
        const vertices = [];
        const faces = [];

        // Tip of the beak
        // Vertex: x, y, z, r, g, b, u, v, nx, ny, nz
        vertices.push(0, 0, 0, color[0], color[1], color[2], 0, 0, 0, 0, 1);

        // Build the beak with semi-circular cross-sections
        for (let i = 1; i <= segments; i++) {
            // Key functions (remain the same)
            const t = i / segments;
            const radiusScale = Math.sqrt(t);
            const currentZ = -length * t;

            // Loop from 0 to segments (inclusive) to create segments+1 vertices
            for (let j = 0; j <= segments; j++) {
                const theta = (j / segments) * Math.PI;

                const x = width * radiusScale * Math.cos(theta);
                const y = thickness * radiusScale * Math.sin(theta); // This will only be >= 0

                let nx = 0, ny = 0, nz = 0;

                if (j === 0 || j === segments) {
                    // This is a vertex on the flat "cap" (the mouth line)
                    // Normal points straight down.
                    nx = 0;
                    ny = -1;
                    nz = 0;
                } else {
                    // This is a vertex on the curved "top" of the beak
                    // Normal points radially outward in XY, similar to generateBeak
                    nx = x;
                    ny = y;
                    nz = 0;
                    const mag = Math.sqrt(nx * nx + ny * ny);
                    if (mag > 0) {
                        nx /= mag;
                        ny /= mag;
                    }
                }

                // Push all 11 vertex components
                vertices.push(x, y, currentZ, color[0], color[1], color[2], t, j / segments, nx, ny, nz);
            }
        }

        // --- Face generation is unchanged ---

        // Create faces for the curved tip
        for (let j = 1; j <= segments; j++) {
            faces.push(0, j, j + 1);
        }

        // Create faces for the curved sides
        for (let i = 0; i < segments - 1; i++) {
            const ring1_start = 1 + i * (segments + 1);
            const ring2_start = 1 + (i + 1) * (segments + 1);

            for (let j = 0; j < segments; j++) {
                const p1 = ring1_start + j;
                const p2 = ring1_start + j + 1;
                const p3 = ring2_start + j;
                const p4 = ring2_start + j + 1;
                faces.push(p1, p3, p2, p2, p3, p4);
            }
        }

        // Create faces for the flat "bottom" cap (on the Y=0 plane)
        for (let i = 0; i < segments; i++) {
            const ring1_start = 1 + (i-1) * (segments + 1);
            const ring2_start = 1 + i * (segments + 1);

            const p1 = (i === 0) ? 0 : ring1_start;
            const p2 = (i === 0) ? 0 : ring1_start + segments;
            const p3 = ring2_start;
            const p4 = ring2_start + segments;

            if (i === 0) {
                // Create the first triangle at the tip
                // Winding order (p1, p4, p3) makes it face down (0, -1, 0)
                faces.push(p1, p4, p3);
            } else {
                // Create the quads for the rest of the cap
                // Winding order (p1, p2, p3) and (p2, p4, p3)
                faces.push(p1, p2, p3, p2, p4, p3);
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

    generateCircle: function(radius, segments, color) {
        const vertices = [];
        const faces = [];
        // Add the center vertex at (0, 0, 0)
        // Vertex format: x, y, z, r, g, b, u, v, nx, ny, nz
        // Normal points up (0, 0, 1) assuming circle is in XY plane
        vertices.push(0, 0, 0, color[0], color[1], color[2], 0.5, 0.5, 0, 0, 1);

        // Add vertices for the circumference
        for (let i = 0; i <= segments; i++) {
            const theta = (i / segments) * 2 * Math.PI;
            const x = radius * Math.cos(theta);
            const y = radius * Math.sin(theta);

            // Texture coordinates (map circular coords to square UV space)
            const u = (x / radius + 1) / 2;
            const v = (y / radius + 1) / 2;

            vertices.push(x, y, 0, color[0], color[1], color[2], u, v, 0, 0, 1);
        }
        // Create the faces using the center vertex (index 0)
        for (let i = 1; i <= segments; i++) {
            const centerIndex = 0;
            const currentIndex = i;
            const nextIndex = i + 1;
            faces.push(centerIndex, currentIndex, nextIndex);
        }
        return { vertices, faces };
    },

    generateAngryEye: function(a, b, c, stack, step, verticalSweepDegrees, color) {
        const vertices = [];
        const faces = [];

        const verticalSweepRad = (verticalSweepDegrees * Math.PI) / 180;

        // The sweep starts from the bottom pole (theta=PI) and goes upwards.
        const endTheta = Math.PI - verticalSweepRad;

        for (let i = 0; i <= stack; i++) {
            for (let j = 0; j <= step; j++) {
                const u = i / stack; // progression along the sweep
                const v = j / step; // progression around the circle

                const theta = Math.PI - u * verticalSweepRad;
                const phi = v * 2 * Math.PI;

                const x = a * Math.sin(theta) * Math.cos(phi);
                const y = b * Math.cos(theta);
                const z = c * Math.sin(theta) * Math.sin(phi);

                // Normal for a partial sphere is its position vector
                const normal = [x, y, z];

                vertices.push(x, y, z, color[0], color[1], color[2], v, u, normal[0], normal[1], normal[2]);
            }
        }

        // Create faces for the curved surface
        for (let i = 0; i < stack; i++) {
            for (let j = 0; j < step; j++) {
                const p1 = i * (step + 1) + j;
                const p2 = p1 + 1;
                const p3 = p1 + (step + 1);
                const p4 = p3 + 1;
                faces.push(p1, p2, p4, p1, p4, p3);
            }
        }

        // If the sweep is less than 180, create a flat top cap to close the shape
        if (verticalSweepDegrees < 180 && verticalSweepDegrees > 0) {
            const capCenterY = b * Math.cos(endTheta);
            const centerIndex = vertices.length / 11; // index for the new center vertex (11 floats)

            // Normal for the cap points straight up (in local Y)
            vertices.push(0, capCenterY, 0, color[0], color[1], color[2], 0.5, 0.5, 0, 1, 0);

            // Create the fan faces for the cap
            const lastRingStartIndex = stack * (step + 1);
            for (let j = 0; j < step; j++) {
                const p1 = lastRingStartIndex + j;
                const p2 = lastRingStartIndex + j + 1;
                // The order (p1, p2, centerIndex) should make the face point upwards
                faces.push(p1, p2, centerIndex);
            }
        }

        return { vertices, faces };
    },

    generateCone: function(radius, height, segments, color) {
        const vertices = [];
        const faces = [];
        const halfHeight = height / 2;

        // --- Vertices ---
        // 1. Tip vertex
        // [x, y, z, r, g, b, u, v, nx, ny, nz]
        vertices.push(0, halfHeight, 0, color[0], color[1], color[2], 0.5, 1, 0, 1, 0); // Normal up

        // 2. Base center vertex (for the bottom cap)
        vertices.push(0, -halfHeight, 0, color[0], color[1], color[2], 0.5, 0, 0, -1, 0); // Normal down

        // 3. Base circumference vertices
        for (let i = 0; i <= segments; i++) {
            const theta = (i / segments) * 2 * Math.PI;
            const x = radius * Math.cos(theta);
            const z = radius * Math.sin(theta);

            // UV coordinates mapping the circle to a square
            const u = (x / radius + 1) / 2;
            const v = (z / radius + 1) / 2;

            // Normal for the side of the cone
            let nx = x;
            let ny = radius / height; // component from the slope
            let nz = z;
            let mag = Math.sqrt(nx*nx + ny*ny + nz*nz);
            if (mag > 0) { nx /= mag; ny /= mag; nz /= mag; }

            vertices.push(x, -halfHeight, z, color[0], color[1], color[2], u, v, nx, ny, nz);
        }

        // --- Faces ---
        const tipIndex = 0;
        const baseCenterIndex = 1;

        for (let i = 0; i < segments; i++) {
            const currentIndex = i + 2;
            const nextIndex = i + 3;

            // Side face
            faces.push(tipIndex, currentIndex, nextIndex);

            // Base face (note the winding order is reversed to face down)
            faces.push(baseCenterIndex, nextIndex, currentIndex);
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

    generateIrregularExtrudedPolygon: function(sides, radius, height, color, irregularityFactor = 0.2) {
        const vertices = [];
        const faces = [];
        const numVerticesPerRing = sides; // Jumlah titik per cincin (atas/bawah)

        // 1. Buat Poin-poin Dasar Poligon (Top Ring)
        const topRingPoints = [];
        for (let i = 0; i < sides; i++) {
            const angle = (i / sides) * 2 * Math.PI;

            // Tambahkan irreguleritas (acak) pada radius setiap titik
            const currentRadius = radius + (Math.random() - 0.5) * radius * irregularityFactor;

            const x = currentRadius * Math.cos(angle);
            const z = currentRadius * Math.sin(angle);
            topRingPoints.push([x, z]);
        }

        // 2. Vertices untuk Bagian ATAS (Top Cap)
        // Center vertex atas
        vertices.push(0, height / 2, 0, color[0], color[1], color[2], 0.5, 0.5, 0, 1, 0); // Normal ke atas (Y+)
        const topCenterIndex = 0;

        // Vertices di pinggir atas
        for (let i = 0; i < sides; i++) {
            const p = topRingPoints[i];
            // Normal untuk top cap adalah (0,1,0)
            vertices.push(p[0], height / 2, p[1], color[0], color[1], color[2], (p[0]/radius+1)/2, (p[1]/radius+1)/2, 0, 1, 0);
        }

        // 3. Faces untuk Bagian ATAS (Top Cap)
        for (let i = 0; i < sides; i++) {
            const currentEdgeIndex = topCenterIndex + 1 + i;
            const nextEdgeIndex = topCenterIndex + 1 + ((i + 1) % sides);
            faces.push(topCenterIndex, currentEdgeIndex, nextEdgeIndex);
        }

        // 4. Vertices untuk Bagian BAWAH (Bottom Cap)
        // Center vertex bawah
        vertices.push(0, -height / 2, 0, color[0], color[1], color[2], 0.5, 0.5, 0, -1, 0); // Normal ke bawah (Y-)
        const bottomCenterIndex = vertices.length / 11 - 1; // Index setelah semua vertex atas

        // Vertices di pinggir bawah (mirip dengan atas, tapi Y-nya negatif)
        for (let i = 0; i < sides; i++) {
            const p = topRingPoints[i]; // Gunakan poin yang sama agar bentuknya lurus
            // Normal untuk bottom cap adalah (0,-1,0)
            vertices.push(p[0], -height / 2, p[1], color[0], color[1], color[2], (p[0]/radius+1)/2, (p[1]/radius+1)/2, 0, -1, 0);
        }

        // 5. Faces untuk Bagian BAWAH (Bottom Cap)
        for (let i = 0; i < sides; i++) {
            const currentEdgeIndex = bottomCenterIndex + 1 + i;
            const nextEdgeIndex = bottomCenterIndex + 1 + ((i + 1) % sides);
            // Urutan harus terbalik agar normal menghadap ke bawah
            faces.push(bottomCenterIndex, nextEdgeIndex, currentEdgeIndex);
        }

        // 6. Vertices dan Faces untuk SISI SAMPING (Side Walls)
        for (let i = 0; i < sides; i++) {
            const p = topRingPoints[i];
            const p_next = topRingPoints[(i + 1) % sides];

            const topVertexIndex = topCenterIndex + 1 + i;
            const bottomVertexIndex = bottomCenterIndex + 1 + i;

            const topNextVertexIndex = topCenterIndex + 1 + ((i + 1) % sides);
            const bottomNextVertexIndex = bottomCenterIndex + 1 + ((i + 1) % sides);

            // Hitung normal untuk sisi samping (vector tegak lurus ke garis pinggir)
            // Ini akan sedikit tricky karena titiknya tidak simetris.
            // Kita bisa ambil normalnya tegak lurus ke segmen garis p - p_next
            const edgeVecX = p_next[0] - p[0];
            const edgeVecZ = p_next[1] - p[1];

            // Normal di XY plane tegak lurus dengan edgeVec adalah (-edgeVecZ, edgeVecX)
            let normalX = -edgeVecZ;
            let normalZ = edgeVecX;
            let normalY = 0; // Karena ini dinding vertikal

            const magN = Math.sqrt(normalX*normalX + normalY*normalY + normalZ*normalZ);
            if (magN > 0) {
                normalX /= magN;
                normalZ /= magN;
            }

            // Normal untuk setiap vertex di sisi samping
            // Untuk memastikan konsistensi, kita perlu menghitung rata-rata normal dari face yang berbagi vertex
            // Namun, untuk sederhana, kita bisa gunakan normal dari segmen garis saat ini untuk 4 vertex ini
            // Ini mungkin sedikit kurang halus di sudut, tapi cukup untuk demonstrasi
            const startIndexSideVerts = vertices.length / 11;

            // Vertex 1 (Top, Current)
            vertices.push(p[0], height / 2, p[1], color[0], color[1], color[2], 0, 1, normalX, normalY, normalZ);
            // Vertex 2 (Top, Next)
            vertices.push(p_next[0], height / 2, p_next[1], color[0], color[1], color[2], 1, 1, normalX, normalY, normalZ);
            // Vertex 3 (Bottom, Current)
            vertices.push(p[0], -height / 2, p[1], color[0], color[1], color[2], 0, 0, normalX, normalY, normalZ);
            // Vertex 4 (Bottom, Next)
            vertices.push(p_next[0], -height / 2, p_next[1], color[0], color[1], color[2], 1, 0, normalX, normalY, normalZ);

            // Faces untuk satu sisi samping (quad)
            // (v1, v3, v4) (v1, v4, v2)
            // Perhatikan bahwa kita tidak menggunakan index dari cap untuk sisi. Kita membuat vertex baru agar normalnya terpisah.
            faces.push(startIndexSideVerts + 0, startIndexSideVerts + 2, startIndexSideVerts + 3); // Triangle 1
            faces.push(startIndexSideVerts + 0, startIndexSideVerts + 3, startIndexSideVerts + 1); // Triangle 2
        }


        return { vertices, faces };
    },

    // --- Fungsi Baru: Bidang Air ---
    generateWaterPlane: function(width, depth, segmentsW, segmentsD, waterColor) {
        const vertices = [];
        const faces = [];
        const halfWidth = width / 2;
        const halfDepth = depth / 2;
        const segmentWidth = width / segmentsW;
        const segmentDepth = depth / segmentsD;

        for (let i = 0; i <= segmentsW; i++) {
            for (let j = 0; j <= segmentsD; j++) {
                const x = i * segmentWidth - halfWidth;
                const z = j * segmentDepth - halfDepth;
                const y = 0; // Bidang datar di Y=0

                // Normal selalu ke atas untuk bidang horizontal (0, 1, 0)
                const nx = 0;
                const ny = 1;
                const nz = 0;

                // UV coordinates
                const u = i / segmentsW;
                const v = j / segmentsD;

                // x, y, z, r, g, b, u, v, nx, ny, nz
                vertices.push(x, y, z, waterColor[0], waterColor[1], waterColor[2], u, v, nx, ny, nz);
            }
        }

        // Generate Faces
        for (let i = 0; i < segmentsW; i++) {
            for (let j = 0; j < segmentsD; j++) {
                const p1 = i * (segmentsD + 1) + j;
                const p2 = p1 + 1;
                const p3 = (i + 1) * (segmentsD + 1) + j;
                const p4 = p3 + 1;

                faces.push(p1, p2, p4);
                faces.push(p1, p4, p3);
            }
        }

        return { vertices, faces };
    },

    generateSkyboxCube: function(size) {
        const s = size / 2;
        const vertices = [
            -s,  s, -s, -s, -s, -s,  s, -s, -s,  s, -s, -s,  s,  s, -s, -s,  s, -s, // Front
            -s, -s,  s, -s,  s,  s,  s,  s,  s,  s,  s,  s,  s, -s,  s, -s, -s,  s, // Back
            -s,  s, -s, -s,  s,  s, -s, -s,  s, -s, -s,  s, -s, -s, -s, -s,  s, -s, // Left
            s, -s, -s,  s, -s,  s,  s,  s,  s,  s,  s,  s,  s,  s, -s,  s, -s, -s, // Right
            -s, -s, -s,  s, -s, -s,  s, -s,  s,  s, -s,  s, -s, -s,  s, -s, -s, -s, // Bottom
            -s,  s, -s, -s,  s,  s,  s,  s,  s,  s,  s,  s,  s,  s, -s, -s,  s, -s  // Top
        ];

        const faces = [
            0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11,
            12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
            24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35
        ];

        return { vertices, faces };
    },

    generateQuad: function(width, height, color) {
        const vertices = [];
        const faces = [];
        const halfWidth = width / 2;
        const halfHeight = height / 2;

        // Define the 4 vertices of the quad in the XY plane
        // [x, y, z, r, g, b, u, v, nx, ny, nz]
        // Normal points forward (+Z)
        vertices.push(-halfWidth, -halfHeight, 0, color[0], color[1], color[2], 0, 0, 0, 0, 1); // Bottom-left
        vertices.push( halfWidth, -halfHeight, 0, color[0], color[1], color[2], 1, 0, 0, 0, 1); // Bottom-right
        vertices.push( halfWidth,  halfHeight, 0, color[0], color[1], color[2], 1, 1, 0, 0, 1); // Top-right
        vertices.push(-halfWidth,  halfHeight, 0, color[0], color[1], color[2], 0, 1, 0, 0, 1); // Top-left

        // Define the 2 triangles that make the quad
        faces.push(0, 1, 2); // Triangle 1
        faces.push(0, 2, 3); // Triangle 2

        return { vertices, faces };
    },
};

export const Geometry = Geometries;