// --- UTILITY FOR GEOMETRY ---
// We place the geometry generation logic in its own object to keep things organized.
(function() {

// (Geometry functions: generateSphere, generateBeak, etc. are all unchanged)
// ... (Your entire Geometry object goes here, unchanged. I've omitted it for brevity.)
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
    };


// --- NEW ModelNode CLASS ---
// This class replaces the old PrinplupPart class
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

        // Function to add a child node
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

        // Set this node's *base* transform (e.g., the flipper's position relative to the body)
        setBaseTransform(matrix) {
            this.baseMatrix = matrix;
            this.localMatrix = matrix; // By default, local = base
        }

        // Set this node's *animated* local transform
        setLocalTransform(matrix) {
            this.localMatrix = matrix;
        }

        // Recursive function to update all matrices in the tree
        updateWorldMatrix(parentWorldMatrix) {
            if (parentWorldMatrix) {
                this.worldMatrix = LIBS.multiply(this.localMatrix, parentWorldMatrix);
            } else {
                this.worldMatrix = this.localMatrix;
            }

            // Now, recursively update all children
            for (const child of this.children) {
                child.updateWorldMatrix(this.worldMatrix);
            }
        }

        // Recursive function to draw this node and all its children
        draw(shader) {
            // Draw ourself (if we have geometry)
            if (this.buffers) {
                const gl = this.gl;
                gl.uniformMatrix4fv(shader.locations.Mmatrix, false, this.worldMatrix);

                gl.bindBuffer(gl.ARRAY_BUFFER, this.buffers.vertex);

                // Stride is now 11 floats (3 pos, 3 color, 2 uv, 3 normal)
                const stride = 4 * (3 + 3 + 2 + 3);
                gl.vertexAttribPointer(shader.locations.position, 3, gl.FLOAT, false, stride, 0);
                gl.vertexAttribPointer(shader.locations.color, 3, gl.FLOAT, false, stride, 3 * 4);
                gl.vertexAttribPointer(shader.locations.texcoord, 2, gl.FLOAT, false, stride, 6 * 4);
                gl.vertexAttribPointer(shader.locations.normal, 3, gl.FLOAT, false, stride, 8 * 4);


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

            // Now, recursively draw all children
            for (const child of this.children) {
                child.draw(shader);
            }
        }
    }


    class Piplup {
        constructor(gl, renderer) {
            this.gl = gl;
            this.renderer = renderer;

            // This is the single root of our scene graph
            this.rootNode = new ModelNode(gl);
            this.modelMatrix = LIBS.get_I4(); // This matrix will control the entire Piplup's rotation

            // MODIFIED: Store references to all animated parts
            this.animatedNodes = {
                bodyNode: null,     // NEW: This node will handle the Y-translation (bob)
                bodyGeometry: null,   // NEW: This node will handle the scale (squash)
                head: null,
                topBeak: null,
                bottomBeak: null,
                leftFlipper: null,
                rightFlipper: null,
                leftLeg: null,
                rightLeg: null
            };

            // NEW: Store the base (default) transformations for animated parts
            this.baseTransforms = {
                topBeak: LIBS.get_I4(),    // NEW
                bottomBeak: LIBS.get_I4(),  // NEW
                leftFlipper: LIBS.get_I4(),
                rightFlipper: LIBS.get_I4(),
                leftLeg: LIBS.get_I4(),
                rightLeg: LIBS.get_I4(),
            };

            this.initParts();
        }

        initParts() {
            const gl = this.gl;
            // Piplup Colors
            const C = {
                BODY: [0.52, 0.80, 1.00], HEAD: [0.20, 0.38, 0.64], BEAK: [1.00, 0.84, 0.00], CAPE: [0.24, 0.42, 0.96],
                EYE_W: [1.00, 1.00, 1.00], BLACK: [0.00, 0.00, 0.00], FEET: [1.00, 0.65, 0.00], WHITE: [1.00, 1.00, 1.00]
            };
            // Helper function to create a translation matrix using your libs.js functions
            const createTransform = (x, y, z) => {
                const m = LIBS.get_I4();
                LIBS.translateX(m, x);
                LIBS.translateY(m, y);
                LIBS.translateZ(m, z);
                return m;
            };

            const createOrientedTransform = (radX, radY, radZ, x, y) => {
                // A. Calculate the precise 'z' on the surface
                const termX = (x * x) / (radX * radX);
                const termY = (y * y) / (radY * radY);
                const z_on_surface = radZ * Math.sqrt(1.0 - termX - termY);

                // Add a tiny offset to prevent clipping
                const z_final = z_on_surface + 0.01;

                // B. Calculate rotation angles to match the surface curve
                const angleY = Math.atan2(x, z_final);
                const angleX = -Math.atan2(y, z_final);

                // C. Build the transformation matrix
                const m = LIBS.get_I4(); // Start with a fresh identity matrix

                // Apply rotations (order matters: Y-axis first, then X-axis)
                LIBS.rotateY(m, angleY);
                LIBS.rotateX(m, angleX);

                // Apply the final position
                LIBS.set_position(m, x, y, z_final);

                return m; // Return the finished matrix
            };

            const headTexture = this.renderer.loadTexture("Resource/piplup_head_texture.png");

            // --- Build the Hierarchy ---
            // All transforms are now RELATIVE to their parent.

// 1. Body (Child of the root)
            // MODIFIED: Create an empty parent node for bobbing
            const bodyNode = new ModelNode(gl);
            bodyNode.setLocalTransform(LIBS.get_I4()); // This node will be animated
            this.rootNode.addChild(bodyNode);
            this.animatedNodes.bodyNode = bodyNode; // Save for animation

            // NEW: Create the visible body geometry as a child of the bobber
            const bodyGeometry = new ModelNode(gl, Geometry.generateSphere(0.8, 0.9, 0.8, 20, 20, C.BODY));
            bodyGeometry.setLocalTransform(LIBS.get_I4()); // This node will be scaled
            bodyNode.addChild(bodyGeometry); // Add to the bobber
            this.animatedNodes.bodyGeometry = bodyGeometry; // Save for animation

            // --- Body Decorations (Children of the Body) ---
            // (Assuming 'bodyNode' is your main body ModelNode)

            // 1. Left White Circle
            const leftDecoCircle = new ModelNode(gl, Geometry.generateCircle(0.2, 20, C.WHITE));
            // We use the original function to get the correct transformation matrix
            const leftDecoTransform = createOrientedTransform(0.8, 1.1, 0.78, -0.3, 0.25);
            leftDecoCircle.setLocalTransform(leftDecoTransform);
            bodyGeometry.addChild(leftDecoCircle);

            // 2. Right White Circle
            const rightDecoCircle = new ModelNode(gl, Geometry.generateCircle(0.2, 20, C.WHITE));
            const rightDecoTransform = createOrientedTransform(0.8, 1.1, 0.78, 0.3, 0.25);
            rightDecoCircle.setLocalTransform(rightDecoTransform);
            bodyGeometry.addChild(rightDecoCircle);

            // 2. Head (Child of the Body)
            const headNode = new ModelNode(gl, Geometry.generateSphere(0.8, 0.8, 0.8, 20, 20, C.HEAD, 0.25), headTexture);
            const headTransform = createTransform(0, 1.5, 0); // Head is 1.5 units *above the body*
            headNode.setLocalTransform(headTransform);
            bodyNode.addChild(headNode);
            this.animatedNodes.head = headNode;

            // 3. Eyes (Children of the Head)
            // Eye position is now relative to the Head (at 0, 1.5, 0)
            // Original eye: (-0.5, 1.4, 0.60)
            // New relative eye: (-0.5, 1.4 - 1.5, 0.60) -> (-0.5, -0.1, 0.60)
            const leftEyeNode = new ModelNode(gl, Geometry.generateSphere(0.1, 0.2, 0.1, 10, 10, C.BLACK));
            leftEyeNode.setLocalTransform(createTransform(-0.5, -0.1, 0.60));
            headNode.addChild(leftEyeNode);

            const leftEyeLightNode = new ModelNode(gl, Geometry.generateSphere(0.05, 0.05, 0.05, 10, 10, C.WHITE));
            leftEyeLightNode.setLocalTransform(createTransform(-0.5, 0, 0.65));
            headNode.addChild(leftEyeLightNode);

            // Right eye
            const rightEyeNode = new ModelNode(gl, Geometry.generateSphere(0.1, 0.2, 0.1, 10, 10, C.BLACK));
            rightEyeNode.setLocalTransform(createTransform(0.5, -0.1, 0.60));
            headNode.addChild(rightEyeNode);

            const rightEyeLightNode = new ModelNode(gl, Geometry.generateSphere(0.05, 0.05, 0.05, 10, 10, C.WHITE));
            rightEyeLightNode.setLocalTransform(createTransform(0.5, 0, 0.65));
            headNode.addChild(rightEyeLightNode);


            // 4. Beak (Child of the Head)
            // Original beak: (0, 1.2, 1.2)
            // New relative beak: (0, 1.2 - 1.5, 1.2) -> (0, -0.3, 1.2)
            const beakUpNode = new ModelNode(gl, Geometry.generateBeakHalf(0.25, 0.3, 0.3, 15, C.BEAK));
            let beakUpMatrix = createTransform(0, -0.35, 1);
            LIBS.rotateX(beakUpMatrix, LIBS.degToRad(15));
            beakUpNode.setLocalTransform(beakUpMatrix);
            headNode.addChild(beakUpNode);
            this.animatedNodes.topBeak = beakUpNode;
            this.baseTransforms.topBeak = beakUpMatrix;

            const beakDownNode = new ModelNode(gl, Geometry.generateBeakHalf(0.25, 0.25, 0.3, 15, C.BEAK));
            let beakDownMatrix = createTransform(0, -0.35, 1);
            LIBS.rotateX(beakDownMatrix, LIBS.degToRad(-25));
            LIBS.rotateZ(beakDownMatrix, LIBS.degToRad(180));
            beakDownNode.setLocalTransform(beakDownMatrix);
            headNode.addChild(beakDownNode);
            this.animatedNodes.bottomBeak = beakDownNode;
            this.baseTransforms.bottomBeak = beakDownMatrix;

            // 5. Hands (Children of the Body)
            const leftHandNode = new ModelNode(gl, Geometry.generateSphere(0.2, 0.7, 0.5, 15, 15, C.BODY));
            let leftHandMatrix = createTransform(-0.8, 0.1, 0.1);
            LIBS.rotateZ(leftHandMatrix, LIBS.degToRad(-20));
            LIBS.rotateX(leftHandMatrix, LIBS.degToRad(-10));
            leftHandNode.setLocalTransform(leftHandMatrix);
            this.baseTransforms.leftFlipper = leftHandMatrix; // NEW: Store base pose
            bodyNode.addChild(leftHandNode);
            this.animatedNodes.leftFlipper = leftHandNode; // Save for animation

            const rightHandNode = new ModelNode(gl, Geometry.generateSphere(0.2, 0.7, 0.5, 15, 15, C.BODY));
            let rightHandMatrix = createTransform(0.8, 0.1, 0.1);
            LIBS.rotateZ(rightHandMatrix, LIBS.degToRad(20));
            LIBS.rotateX(rightHandMatrix, LIBS.degToRad(-10));
            rightHandNode.setLocalTransform(rightHandMatrix);
            this.baseTransforms.rightFlipper = rightHandMatrix; // NEW: Store base pose
            bodyNode.addChild(rightHandNode);
            this.animatedNodes.rightFlipper = rightHandNode; // MODIFIED: Corrected typo (was leftFlipper)

            // ... (add right flipper as a child of bodyNode) ...

            // --- Cape and Back Details (Children of the Body) ---
            // (Assuming 'bodyNode' is the ModelNode for the Piplup's body)

            // --- Legs and Feet (Children of Body) ---

            // (Assuming 'bodyNode' is your main body ModelNode)
            // (And C.FEET and C.BODY are defined in your colors)

            // 1. Left Leg
            const leftLegNode = new ModelNode(gl, Geometry.generateSphere(0.25, 0.5, 0.25, 10, 10, C.BODY));
            const leftLegTransform = createTransform(-0.4, -0.5, -0.1); // NEW: Capture transform
            leftLegNode.setLocalTransform(leftLegTransform);
            this.baseTransforms.leftLeg = leftLegTransform; // NEW: Store base pose
            bodyNode.addChild(leftLegNode);
            this.animatedNodes.leftLeg = leftLegNode; // NEW: Save for animation


            // 2. Left Foot (as a child of the Left Leg)
            const leftFootNode = new ModelNode(gl, Geometry.generateSphere(0.25, 0.15, 0.45, 10, 10, C.FEET));
            // The transform is (0, -0.5, 0.2) *relative* to the leg's position
            leftFootNode.setLocalTransform(createTransform(0, -0.5, 0.3));
            leftLegNode.addChild(leftFootNode); // <-- Add to leg, NOT body

            // 3. Right Leg
            const rightLegNode = new ModelNode(gl, Geometry.generateSphere(0.25, 0.5, 0.25, 10, 10, C.BODY));
            const rightLegTransform = createTransform(0.4, -0.5, -0.1); // NEW: Capture transform
            rightLegNode.setLocalTransform(rightLegTransform);
            this.baseTransforms.rightLeg = rightLegTransform; // NEW: Store base pose
            bodyNode.addChild(rightLegNode);
            this.animatedNodes.rightLeg = rightLegNode; // NEW: Save for animation

            // 4. Right Foot (as a child of the Right Leg)
            const rightFootNode = new ModelNode(gl, Geometry.generateSphere(0.25, 0.15, 0.45, 10, 10, C.FEET));
            // The transform is (0, -0.5, 0.2) *relative* to the leg's position
            rightFootNode.setLocalTransform(createTransform(0, -0.5, 0.3));
            rightLegNode.addChild(rightFootNode); // <-- Add to leg, NOT body

            // 1. The Spline Tube (Main Cape/Tail)
            const capeSplineNode = new ModelNode(gl, Geometry.generateTubeFromSpline(
                [
                    [0.0, 0.0, -0.5],   // Start point on the lower back
                    [0.0, 0.05, -1.5],
                    [0.0, -0.2, -1.3], // First curve point
                    [0.0, -0.4, -1.8],  // Second curve point
                    [0.0, -0.7, -1.5]   // End point, slightly flared out
                ],
                100, // Segments
                0.15, // Radius
                20,   // Radial segments
                [0.20, 0.38, 0.64] // Color
            ));
            capeSplineNode.setLocalTransform(createTransform(0, 0, 0.5));
            bodyNode.addChild(capeSplineNode);

            // 2. Right Shoulder Sphere
            const capeRightShoulder = new ModelNode(gl, Geometry.generateSphere(0.5, 0.2, 0.4, 20, 20, C.HEAD));
            let m_rightShoulder = createTransform(0.3, 0.8, 0.4);
            LIBS.rotateY(m_rightShoulder, LIBS.degToRad(-60));
            LIBS.rotateZ(m_rightShoulder, LIBS.degToRad(-30));
            capeRightShoulder.setLocalTransform(m_rightShoulder);
            bodyNode.addChild(capeRightShoulder);

            // 3. Left Shoulder Sphere
            const capeLeftShoulder = new ModelNode(gl, Geometry.generateSphere(0.5, 0.2, 0.4, 20, 20, C.HEAD));
            let m_leftShoulder = createTransform(-0.3, 0.8, 0.4);
            LIBS.rotateY(m_leftShoulder, LIBS.degToRad(60));
            LIBS.rotateZ(m_leftShoulder, LIBS.degToRad(30));
            capeLeftShoulder.setLocalTransform(m_leftShoulder);
            bodyNode.addChild(capeLeftShoulder);

            // 4. Neck Ring Sphere
            const capeNeckRing = new ModelNode(gl, Geometry.generateSphere(0.7, 0.12, 0.6, 20, 20, C.HEAD));
            let m_neckRing = createTransform(0, 0.8, 0);
            // (Original rotations were commented out, so I've left them out)
            capeNeckRing.setLocalTransform(m_neckRing);
            bodyNode.addChild(capeNeckRing);

            // 5. Back Sphere
            const capeBack = new ModelNode(gl, Geometry.generateSphere(1, 1.1, 0.25, 20, 20, C.HEAD));
            let m_back = createTransform(0, 0, -0.7);
            LIBS.rotateX(m_back, LIBS.degToRad(20));
            capeBack.setLocalTransform(m_back);
            bodyNode.addChild(capeBack);

            // ... (Continue for all other parts, attaching them to either the body or head) ...

        }

        // NEW: Function to update animations
        updateAnimation(time) {
            const speed = 1; // Controls the speed of the run
            const bobAmount = Math.abs(Math.sin(time * speed) * 0.5) * 0.2 - 0.05; // Bob up and down
            const swingAngle = Math.sin(time * speed) * 0.3 - 0.4; // Swing angle in radians (~40 deg)
            const legSwingAngle = Math.sin(time * speed) * 0.5; // Legs swing a bit less

            // 1. Body Bob (with Squash and Stretch)

            // --- Calculate Squash ---
            const squashFactor = Math.abs(Math.sin(time * speed) * 0.5);
            const squashAmount = 0.1;
            const scaleY = 1.0 - squashFactor * squashAmount;
            const scaleXZ = 1.0 + squashFactor * squashAmount * 0.5;

            // --- Create Translation Matrix ---
            const T_body = LIBS.get_I4();
            LIBS.translateY(T_body, bobAmount + 1.6);

            // --- Create Scale Matrix ---
            const S_body = LIBS.get_I4();
            LIBS.scaleY(S_body, scaleY);
            LIBS.scaleX(S_body, scaleXZ);
            LIBS.scaleZ(S_body, scaleXZ);

            // --- Apply Matrices ---
            // MODIFIED: Apply Translation to the parent bobber
            this.animatedNodes.bodyNode.setLocalTransform(T_body);

            // MODIFIED: Apply Scale ONLY to the body geometry
            this.animatedNodes.bodyGeometry.setLocalTransform(S_body);


            // 2. Left Flipper (Flaps up with top pivot)
            // Create the animation matrix: M_anim = T(0.7) * R(angle) * T(-0.7)
            const T_flipper_up = LIBS.get_I4();
            LIBS.translateY(T_flipper_up, -0.6);
            LIBS.translateX(T_flipper_up, 0.5);

            const R_flipper = LIBS.get_I4();
            // MODIFIED: Changed from rotateX to rotateZ
            LIBS.rotateZ(R_flipper, swingAngle);

            const T_flipper_down = LIBS.get_I4();
            LIBS.translateY(T_flipper_down, 0.6);
            LIBS.translateX(T_flipper_down, -0.5);

            // Combine them: M_anim = T_up * (R * T_down)
            let leftFlipperAnim = LIBS.multiply(R_flipper, T_flipper_down);
            leftFlipperAnim = LIBS.multiply(T_flipper_up, leftFlipperAnim);

            // Apply animation to the base pose: M_final = M_base * M_anim
            const leftFlipperFinal = LIBS.multiply(this.baseTransforms.leftFlipper, leftFlipperAnim);
            this.animatedNodes.leftFlipper.setLocalTransform(leftFlipperFinal);


            // 3. Right Flipper (Flaps down with top pivot)
// --- NEW: Create mirrored pivot matrices for the right flipper ---
// Pivot point is at (x: 0.5, y: 0.6) relative to its center
            const T_flipper_up_right = LIBS.get_I4();
            LIBS.translateY(T_flipper_up_right, -0.6);
            LIBS.translateX(T_flipper_up_right, -0.5); // Mirrored X

            const R_flipper_right = LIBS.get_I4();
            LIBS.rotateZ(R_flipper_right, swingAngle + 0.8); // Note the negative angle

            const T_flipper_down_right = LIBS.get_I4();
            LIBS.translateY(T_flipper_down_right, 0.6);
            LIBS.translateX(T_flipper_down_right, 0.5); // Mirrored X

// --- MODIFIED: Use the new right-flipper matrices ---
// Combine them: M_anim = T_up_right * (R_right * T_down_right)
            let rightFlipperAnim = LIBS.multiply(R_flipper_right, T_flipper_down_right);
            rightFlipperAnim = LIBS.multiply(T_flipper_up_right, rightFlipperAnim);

// Apply animation to the base pose: M_final = M_base * M_anim
            const rightFlipperFinal = LIBS.multiply(this.baseTransforms.rightFlipper, rightFlipperAnim);
            this.animatedNodes.rightFlipper.setLocalTransform(rightFlipperFinal);

            // 4. Left Leg (Swings backward, opposite of left flipper)
            const leftLegAnim = LIBS.get_I4();
            LIBS.rotateX(leftLegAnim, -legSwingAngle);
            const leftLegFinal = LIBS.multiply(leftLegAnim, this.baseTransforms.leftLeg);
            this.animatedNodes.leftLeg.setLocalTransform(leftLegFinal);

            // 5. Right Leg (Swings forward, opposite of right flipper)
            const rightLegAnim = LIBS.get_I4();
            LIBS.rotateX(rightLegAnim, legSwingAngle);
            const rightLegFinal = LIBS.multiply(rightLegAnim, this.baseTransforms.rightLeg);
            this.animatedNodes.rightLeg.setLocalTransform(rightLegFinal);


// 6. Beak "Talking" Animation
            // Creates a small "chirp" animation.
            // (Math.sin + 1)/2 gives a 0-to-1 range. * 0.2 gives a max 0.2 rad angle.
            const beakAngle = (Math.sin(time * speed * 1.5 + 2) + 1) * 0.5 * 0.2;

            // Define the pivot matrices. Hinge is at z = -0.3
            const T_beak_hinge = LIBS.get_I4();
            LIBS.translateZ(T_beak_hinge, 0.7); // 1. Move hinge to origin

            const T_beak_back = LIBS.get_I4();
            LIBS.translateZ(T_beak_back, -0.7); // 3. Move hinge back

            // --- Top Beak ---
            const R_top_beak = LIBS.get_I4();
            LIBS.rotateX(R_top_beak, -beakAngle); // 2. Rotate (opens up)

            // Combine: M_anim = T_back * R_anim * T_hinge
            let topBeakAnim = LIBS.multiply(R_top_beak, T_beak_hinge);
            topBeakAnim = LIBS.multiply(T_beak_back, topBeakAnim);

            // Apply to base pose: M_final = M_base * M_anim
            const topBeakFinal = LIBS.multiply(this.baseTransforms.topBeak, topBeakAnim);
            this.animatedNodes.topBeak.setLocalTransform(topBeakFinal);

            // --- Bottom Beak ---
            const R_bottom_beak = LIBS.get_I4();
            LIBS.rotateX(R_bottom_beak, beakAngle); // 2. Rotate (opens down)

            // Combine: M_anim = T_back * R_anim * T_hinge
            let bottomBeakAnim = LIBS.multiply(R_bottom_beak, T_beak_hinge);
            bottomBeakAnim = LIBS.multiply(T_beak_back, bottomBeakAnim);

            // Apply to base pose: M_final = M_base * M_anim
            const bottomBeakFinal = LIBS.multiply(this.baseTransforms.bottomBeak, bottomBeakAnim);
            this.animatedNodes.bottomBeak.setLocalTransform(bottomBeakFinal);    }


        draw(shader, parentMatrix) {
            // 1. Update the entire tree's matrices
            // The parentMatrix (from the ice island) is now passed in
            this.rootNode.updateWorldMatrix(parentMatrix);

            // 2. Start the recursive draw
            this.rootNode.draw(shader);
        }}



// --- Prinplup CONTAINER CLASS (Refactored) ---
// Manages all the parts that make up the Prinplup.
    class Prinplup {
        constructor(gl, renderer) {
            this.gl = gl;
            this.renderer = renderer;

            this.rootNode = new ModelNode(gl);

            // MODIFIED: Set the static Y-offset here
            this.modelMatrix = LIBS.get_I4();
            // LIBS.translateY(this.modelMatrix, 2.5); // Prinplup's base offset

            // Store references to nodes that need animation
            this.animatedNodes = {
                breathingNode: null,
                eyeNodes: [],
                leftHand: null,   // NEW
                rightHand: null,   // NEW
            };

            // NEW: Add a place to store base poses
            this.baseTransforms = {
                leftHand: LIBS.get_I4(),
                rightHand: LIBS.get_I4()
            };

            this.initParts();
        }

        initParts() {
            const gl = this.gl;
            const C = {
                BODY: [0.61, 0.84, 0.89], HEAD: [0.20, 0.38, 0.64], BEAK: [1.00, 0.84, 0.00], CAPE: [0.24, 0.42, 0.96],
                EYE_W: [1.00, 1.00, 1.00], BLACK: [0.00, 0.00, 0.00], FEET: [1.00, 0.65, 0.00], WHITE: [1.00, 1.00, 1.00], TAIL: [0.08, 0.32, 0.60]
            };

            const createTransform = (x, y, z) => {
                const m = LIBS.get_I4();
                LIBS.translateX(m, x);
                LIBS.translateY(m, y);
                LIBS.translateZ(m, z);
                return m;
            };

            const createOrientedTransform = (radX, radY, radZ, x, y) => {
                const termX = (x * x) / (radX * radX);
                const termY = (y * y) / (radY * radZ); // Bug? Should be radY * radY
                const z_on_surface = radZ * Math.sqrt(1.0 - termX - termY);
                const z_final = z_on_surface + 0.01;
                const angleY = Math.atan2(x, z_final);
                const angleX = -Math.atan2(y, z_final);
                const m = LIBS.get_I4();
                LIBS.rotateY(m, angleY);
                LIBS.rotateX(m, angleX);
                LIBS.set_position(m, x, y, z_final);
                return m;
            };

            const headTexture = this.renderer.loadTexture("Resource/prinplup_texture_bare.png");
            const handTexture = this.renderer.loadTexture("Resource/prinplup_hand_texture.png");

            const body_profile = [
                [0.0, 1.95, 0], [0.15, 1.95, 0], [0.3, 1.9, 0], [0.4, 1.8, 0], [0.475, 1.7, 0],
                [0.5, 1.6, 0], [0.525, 1.5, 0], [0.53, 1.4, 0], [0.535, 1.3, 0], [0.55, 1.2, 0],
                [0.57, 1.1, 0], [0.59, 1.0, 0], [0.62, 0.9, 0], [0.65, 0.8, 0], [0.70, 0.7, 0],
                [0.75, 0.6, 0], [0.80, 0.5, 0], [0.82, 0.4, 0], [0.84, 0.3, 0], [0.86, 0.2, 0],
                [0.88, 0.1, 0], [0.92, 0, 0], [0.94, -0.1, 0], [0.96, -0.2, 0], [0.96, -0.3, 0],
                [0.96, -0.4, 0], [0.96, -0.5, 0], [0.96, -0.6, 0], [0.96, -0.7, 0], [0.96, -0.8, 0],
                [0.96, -0.9, 0], [0.90, -1.0, 0], [0.84, -1.1, 0], [0.64, -1.2, 0], [0.44, -1.21, 0],
                [0.24, -1.22, 0], [0.04, -1.23, 0], [0.0, -1.24, 0],
            ];

            // --- Build the Hierarchy ---

            // 1. Create the invisible "breathing" node.
            // All parts that move up/down together will be a child of this.
            const breathingNode = new ModelNode(gl);
            this.rootNode.addChild(breathingNode);
            this.animatedNodes.breathingNode = breathingNode;

            // 2. Create and add all parts as children of 'breathingNode'

            // Body
            const bodyNode = new ModelNode(gl, Geometry.generateLathe(body_profile, 30, C.BODY), headTexture);
            const bodyTrans = (() => {
                const m = createTransform(0, 0.2, 0);
                LIBS.rotateY(m, Math.PI / 2);
                LIBS.scale(m, 1.2);
                return m
            })();
            bodyNode.setBaseTransform(bodyTrans);
            breathingNode.addChild(bodyNode);

            // Body Decorations
            const deco1 = new ModelNode(gl, Geometry.generateCircle(0.17, 20, C.WHITE));
            deco1.setBaseTransform(createOrientedTransform(0.6, 1.1, 1.15, -0.4, 0.25));
            breathingNode.addChild(deco1);

            const deco2 = new ModelNode(gl, Geometry.generateCircle(0.17, 20, C.WHITE));
            deco2.setBaseTransform(createOrientedTransform(0.6, 1.1, 1.15, 0.4, 0.25));
            breathingNode.addChild(deco2);

            const deco3 = new ModelNode(gl, Geometry.generateCircle(0.17, 20, C.WHITE));
            deco3.setBaseTransform(createTransform(-0.4, -0.4, 0.94));
            breathingNode.addChild(deco3);

            const deco4 = new ModelNode(gl, Geometry.generateCircle(0.17, 20, C.WHITE));
            deco4.setBaseTransform(createTransform(0.4, -0.4, 0.94));
            breathingNode.addChild(deco4);

            // Head Disks
            const disk1 = new ModelNode(gl, Geometry.generateSphere(0.7, 0.12, 0.55, 20, 20, C.BEAK));
            disk1.setBaseTransform((() => {
                let m = createTransform(-0.2, 2.1, 0);
                LIBS.rotateZ(m, LIBS.degToRad(120));
                return m;
            })());
            breathingNode.addChild(disk1);

            const disk2 = new ModelNode(gl, Geometry.generateSphere(0.7, 0.12, 0.55, 20, 20, C.BEAK));
            disk2.setBaseTransform((() => {
                let m = createTransform(0.2, 2.1, 0);
                LIBS.rotateZ(m, LIBS.degToRad(60));
                return m;
            })());
            breathingNode.addChild(disk2);

            // Eyes (These are animated, so we save them)
            const leftEye1 = new ModelNode(gl, Geometry.generateAngryEye(0.15, 0.2, 0.1, 10, 10, 130, C.WHITE));
            leftEye1.setBaseTransform((() => {
                let m = createTransform(-0.3, 1.9, 0.4);
                LIBS.rotateZ(m, LIBS.degToRad(330));
                return m;
            }) ());
            breathingNode.addChild(leftEye1);
            this.animatedNodes.eyeNodes.push(leftEye1);

            const leftEye2 = new ModelNode(gl, Geometry.generateAngryEye(0.1, 0.15, 0.12, 10, 10, 130, C.HEAD));
            leftEye2.setBaseTransform((() => {
                let m = createTransform(-0.33, 1.9, 0.42);
                LIBS.rotateZ(m, LIBS.degToRad(330));
                return m;
            }) ());
            breathingNode.addChild(leftEye2);
            this.animatedNodes.eyeNodes.push(leftEye2);

            const leftEye3 = new ModelNode(gl, Geometry.generateAngryEye(0.07, 0.12, 0.12, 10, 10, 130, C.BLACK));
            leftEye3.setBaseTransform((() => {
                let m = createTransform(-0.3, 1.9, 0.45);
                LIBS.rotateZ(m, LIBS.degToRad(330));
                return m;
            }) ());
            breathingNode.addChild(leftEye3);
            this.animatedNodes.eyeNodes.push(leftEye3);

            const leftEye4 = new ModelNode(gl, Geometry.generateSphere(0.03, 0.03, 0.07, 10, 10, C.WHITE));
            leftEye4.setBaseTransform(createTransform(-0.3, 1.93, 0.52));
            breathingNode.addChild(leftEye4);
            this.animatedNodes.eyeNodes.push(leftEye4);

            // Right Eye
            const rightEye1 = new ModelNode(gl, Geometry.generateAngryEye(0.15, 0.2, 0.1, 10, 10, 130, C.WHITE));
            rightEye1.setBaseTransform((() => {
                let m = createTransform(0.3, 1.9, 0.4);
                LIBS.rotateZ(m, LIBS.degToRad(30));
                return m;
            }) ());
            breathingNode.addChild(rightEye1);
            this.animatedNodes.eyeNodes.push(rightEye1);

            const rightEye2 = new ModelNode(gl, Geometry.generateAngryEye(0.1, 0.15, 0.12, 10, 10, 130, C.HEAD));
            rightEye2.setBaseTransform((() => {
                let m = createTransform(0.33, 1.9, 0.42);
                LIBS.rotateZ(m, LIBS.degToRad(30));
                return m;
            }) ());
            breathingNode.addChild(rightEye2);
            this.animatedNodes.eyeNodes.push(rightEye2);

            const rightEye3 = new ModelNode(gl, Geometry.generateAngryEye(0.07, 0.12, 0.12, 10, 10, 130, C.BLACK));
            rightEye3.setBaseTransform((() => {
                let m = createTransform(0.3, 1.9, 0.45);
                LIBS.rotateZ(m, LIBS.degToRad(30));
                return m;
            }) ());
            breathingNode.addChild(rightEye3);
            this.animatedNodes.eyeNodes.push(rightEye3);

            const rightEye4 = new ModelNode(gl, Geometry.generateSphere(0.03, 0.03, 0.07, 10, 10, C.WHITE));
            rightEye4.setBaseTransform(createTransform(0.3, 1.93, 0.52));
            breathingNode.addChild(rightEye4);
            this.animatedNodes.eyeNodes.push(rightEye4);

            // Beak
            const beak = new ModelNode(gl, Geometry.generateBeak(0.22, 0.35, 0.55, 20, C.BEAK));
            beak.setBaseTransform((() => {
                let m = createTransform(0, 1.8, 0.9);
                LIBS.rotateX(m, LIBS.degToRad(5));
                return m;
            })());
            breathingNode.addChild(beak);

            const beakCone = new ModelNode(gl, Geometry.generateCone(0.15, 0.3, 15, C.BEAK));
            beakCone.setBaseTransform((() => {
                let m = createTransform(0, 1.95, 0.6);
                LIBS.rotateX(m, LIBS.degToRad(-15));
                return m;
            })());
            breathingNode.addChild(beakCone);

            // Hands
            const leftHand = new ModelNode(gl, Geometry.generateSphere(0.2, 1.5, 0.5, 15, 15, C.TAIL), handTexture);
            // NEW: Store the matrix and the node
            const leftHandMatrix = (() => {
                let m = createTransform(-1.2, 0.42, 0.1);
                LIBS.rotateZ(m, LIBS.degToRad(-40));
                LIBS.rotateX(m, LIBS.degToRad(-10));
                return m;
            })();
            leftHand.setBaseTransform(leftHandMatrix);
            this.baseTransforms.leftHand = leftHandMatrix;
            this.animatedNodes.leftHand = leftHand;
            breathingNode.addChild(leftHand);

            const rightHand = new ModelNode(gl, Geometry.generateSphere(0.2, 1.5, 0.5, 15, 15, C.TAIL), handTexture);
            // NEW: Store the matrix and the node
            const rightHandMatrix = (() => {
                let m = createTransform(1.2, 0.42, 0.1);
                LIBS.rotateZ(m, LIBS.degToRad(40));
                LIBS.rotateX(m, LIBS.degToRad(-10));
                return m;
            })();
            rightHand.setBaseTransform(rightHandMatrix);
            this.baseTransforms.rightHand = rightHandMatrix;
            this.animatedNodes.rightHand = rightHand;
            breathingNode.addChild(rightHand);

            // Legs (Creating parent nodes for feet)
            const leftLeg = new ModelNode(gl, Geometry.generateSphere(0.3, 0.8, 0.4, 10, 10, C.BODY));
            leftLeg.setBaseTransform(createTransform(-0.65, -1.2, 0.1));
            breathingNode.addChild(leftLeg); // Add leg to body

            const leftFoot1 = new ModelNode(gl, Geometry.generateSphere(0.3, 0.12, 0.35, 10, 10, C.FEET));
            leftFoot1.setBaseTransform(createTransform(-0.65, -1.9, 0.24));
            breathingNode.addChild(leftFoot1); // Add foot to body

            const leftFoot2 = new ModelNode(gl, Geometry.generateSphere(0.3, 0.1, 0.5, 10, 10, C.FEET));
            leftFoot2.setBaseTransform(createTransform(-0.65, -2, 0.4));
            breathingNode.addChild(leftFoot2); // Add foot to body

            const rightLeg = new ModelNode(gl, Geometry.generateSphere(0.3, 0.8, 0.4, 10, 10, C.BODY));
            rightLeg.setBaseTransform(createTransform(0.65, -1.2, 0.1));
            breathingNode.addChild(rightLeg);

            const rightFoot1 = new ModelNode(gl, Geometry.generateSphere(0.3, 0.12, 0.35, 10, 10, C.FEET));
            rightFoot1.setBaseTransform(createTransform(0.65, -1.9, 0.24));
            breathingNode.addChild(rightFoot1);

            const rightFoot2 = new ModelNode(gl, Geometry.generateSphere(0.3, 0.1, 0.6, 10, 10, C.FEET));
            rightFoot2.setBaseTransform(createTransform(0.65, -2, 0.4));
            breathingNode.addChild(rightFoot2);

            // Tail
            const tail = new ModelNode(gl, Geometry.generateTaperedShapeFromSpline(
                [[0.0, -0, -0.2], [0.0, -0.9, -1.2], [0.0, -1.2, -2.0]],
                50, [0.9, 0.6], [0.01, 0.01], 20, C.TAIL
            ));
            tail.setBaseTransform(createTransform(0, -0.2, 0.1));
            breathingNode.addChild(tail);
        }

        // NEW: Central animation update function
        updateAnimation(animValues) {
            // 1. Apply 'bodyBreathe' Y-translation to the main breathing node
            const T_breath = LIBS.get_I4();
            LIBS.translateY(T_breath, animValues.body + 2.5);
            this.animatedNodes.breathingNode.setLocalTransform(T_breath);

            // Emang matanya goyang2 ta?
            // // 2. Apply 'eyeBreathe' Y-scale to all eye parts
            // const S_eye = LIBS.get_I4();
            // LIBS.scaleY(S_eye, animValues.eye || 1.0);
            //
            // this.animatedNodes.eyeNodes.forEach(eyeNode => {
            //     // Combine the base pose with the animation: M_local = M_base * M_anim
            //     const finalEyeMatrix = LIBS.multiply(eyeNode.baseMatrix, S_eye);
            //     eyeNode.setLocalTransform(finalEyeMatrix);
            // });

            // (Future animations like 'diskBreathe' would be added here)

            // 3. NEW: Apply 'flapAngle' Z-rotation to hands
            const flapAngle = animValues.flapAngle || 0.0;
            const pivotY = 1.5; // From the hand's geometry radius
            const pivotX = 0.7;

            // Create pivot matrices
            let T_up = LIBS.get_I4();
            LIBS.translateY(T_up, -pivotY);
            LIBS.translateX(T_up, pivotX);
            let T_down = LIBS.get_I4();
            LIBS.translateY(T_down, pivotY);
            LIBS.translateX(T_down, -pivotX);

            // --- Left Hand ---
            const R_left = LIBS.get_I4();
            LIBS.rotateZ(R_left, flapAngle);

            // Combine for animation matrix: M_anim = T_up * R_z * T_down
            let leftAnim = LIBS.multiply(R_left, T_down);
            leftAnim = LIBS.multiply(T_up, leftAnim);

            // Combine with base: M_final = M_base * M_anim
            const leftFinal = LIBS.multiply(this.baseTransforms.leftHand, leftAnim);
            this.animatedNodes.leftHand.setLocalTransform(leftFinal);

            // --- Right Hand ---
            // Create pivot matrices
            T_up = LIBS.get_I4();
            LIBS.translateY(T_up, -pivotY);
            LIBS.translateX(T_up, -pivotX);
            T_down = LIBS.get_I4();
            LIBS.translateY(T_down, pivotY);
            LIBS.translateX(T_down, pivotX);

            const R_right = LIBS.get_I4();
            LIBS.rotateZ(R_right, -flapAngle); // Opposite direction

            // Combine for animation matrix: M_anim = T_up * R_z * T_down
            let rightAnim = LIBS.multiply(R_right, T_down);
            rightAnim = LIBS.multiply(T_up, rightAnim);

            // Combine with base: M_final = M_base * M_anim
            const rightFinal = LIBS.multiply(this.baseTransforms.rightHand, rightAnim);
            this.animatedNodes.rightHand.setLocalTransform(rightFinal);
        }

        draw(shader, parentMatrix) {
            // 'parentMatrix' is the animated matrix from the ice island
            // 'this.modelMatrix' is the mouse-drag rotation
            const finalParentMatrix = LIBS.multiply(parentMatrix, this.modelMatrix);

            // Update all world matrices starting from the root
            this.rootNode.updateWorldMatrix(finalParentMatrix);

            // Start the recursive draw
            this.rootNode.draw(shader);
        }
    }

    class Empoleon {
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
                LIBS.rotateZ(R_left, flapAngle); // Rotate around Z axis

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
                LIBS.rotateZ(R_right, -flapAngle); // Opposite direction

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
// --- ENVIRONMENT CONTAINER CLASS (Refactored) ---
    class Environment {
        constructor(gl, renderer) {
            this.gl = gl;
            this.renderer = renderer;

            this.rootNode = new ModelNode(gl); // Root of the environment scene
            this.modelMatrix = LIBS.get_I4();  // Static matrix for the whole env

            this.animatedNodes = {
                iceIsland: null
            };
            this.animationTime = 0;

            this.initParts();
        }

        initParts() {
            const gl = this.gl;

            const createTransform = (x, y, z) => {
                const m = LIBS.get_I4();
                LIBS.translateX(m, x);
                LIBS.translateY(m, y);
                LIBS.translateZ(m, z);
                return m;
            };

            const C = {
                SNOW_WHITE: [0.95, 0.98, 1.0],
                ICE_BLUE: [0.6, 0.8, 0.95],
                WATER_BLUE: [0.192, 0.502, 0.647]
            };

            // 1. Ice Island
            const iceIslandNode = new ModelNode(gl, Geometry.generateIrregularExtrudedPolygon(8, 15, 1, C.SNOW_WHITE, 0.5));
            iceIslandNode.setBaseTransform(LIBS.get_I4());
            this.rootNode.addChild(iceIslandNode);
            this.animatedNodes.iceIsland = iceIslandNode; // Save for animation

            // 2. Water
            const waterNode = new ModelNode(gl, Geometry.generateWaterPlane(200, 200, 1, 1, C.WATER_BLUE));
            waterNode.setBaseTransform(LIBS.get_I4());
            this.rootNode.addChild(waterNode);
        }

        // NEW: Central animation update
        updateAnimation() {
            this.animationTime += 0.02;

            // 1. Calculate float animation
            const amplitude = 0.1;
            const floatY = Math.sin(this.animationTime) * amplitude;

            // 2. Apply animation to the ice island
            // Note: We combine with the base matrix
            const T_float = LIBS.get_I4();
            LIBS.translateY(T_float, floatY);

            const iceNode = this.animatedNodes.iceIsland;
            const finalIceMatrix = LIBS.multiply(iceNode.baseMatrix, T_float);
            iceNode.setLocalTransform(finalIceMatrix);
        }

        // MODIFIED: Return the final computed world matrix of the island
        getIceIslandWorldMatrix() {
            // Ensure the matrix is up-to-date before returning it
            this.rootNode.updateWorldMatrix(this.modelMatrix);
            return this.animatedNodes.iceIsland.worldMatrix;
        }

        draw(shader) {
            // Update all world matrices
            this.rootNode.updateWorldMatrix(this.modelMatrix);
            // Start recursive draw
            this.rootNode.draw(shader);
        }
    }


// --- RENDERER CLASS (Modified) ---
    class Renderer {
        constructor(canvasId) {
            this.canvas = document.getElementById(canvasId);
            this.canvas.width = window.innerWidth;
            this.canvas.height = window.innerHeight;

            this.gl = this.canvas.getContext("webgl", { antialias: true });
            if (!this.gl) throw new Error("WebGL not supported");

            this.shader = this.createShaderProgram();

            // NEW: Create skybox shader, buffers, and texture
            this.skyboxShader = this.createSkyboxShaderProgram();
            this.skyboxBuffers = this.createSkyboxBuffers();
            this.skyboxTexture = this.loadCubemapTexture([
                'Resource/skybox/frozenrt.png', // Positive X (Right)
                'Resource/skybox/frozenlf.png', // Negative X (Left)
                'Resource/skybox/frozenup.png', // Positive Y (Top)
                'Resource/skybox/frozendn.png', // Negative Y (Bottom)
                'Resource/skybox/frozenbk.png', // Positive Z (Back)
                'Resource/skybox/frozenft.png'  // Negative Z (Front)
            ]);

            // Piplup
            this.piplup = new Piplup(this.gl, this);

            this.Prinplup = new Prinplup(this.gl, this);

            // Empoleon
            this.empoleon = new Empoleon(this.gl, this);

            this.environment = new Environment(this.gl, this);

            // NEW: Create a static matrix for Piplup's offset
            this.piplupModelMatrix = LIBS.get_I4();
            LIBS.translateX(this.piplupModelMatrix, -5.0); // Position 5 units to the left
            LIBS.scale(this.piplupModelMatrix, 0.9); // Make it 90% of the size

            this.viewMatrix = LIBS.get_I4();
            LIBS.translateZ(this.viewMatrix, -20);
            LIBS.translateY(this.viewMatrix, -2);
            this.projMatrix = LIBS.get_projection(40, this.canvas.width / this.canvas.height, 1, 100);

            // Static matrix for Empoleon
            this.empoleonModelMatrix = LIBS.get_I4();
            LIBS.translateX(this.empoleonModelMatrix, 5.0); // Position 5 units to the right
            LIBS.scale(this.empoleonModelMatrix, 1.0); // Keep it at 100% size

            this.initInputHandlers();
            this.startRenderLoop();
        }

        createShaderProgram() {
            const gl = this.gl;
            const vsSource = `
            attribute vec3 position;
            attribute vec3 color;
            attribute vec2 texcoord;
            attribute vec3 normal;
            uniform mat4 Mmatrix, Vmatrix, Pmatrix;
            varying vec3 vColor;
            varying vec2 vTexcoord;
            varying vec3 vNormal;
            varying vec3 v_worldPosition;
            void main(void) {
                gl_Position = Pmatrix * Vmatrix * Mmatrix * vec4(position, 1.);
                vColor = color;
                vTexcoord = texcoord;
                vNormal = mat3(Mmatrix) * normal;
                v_worldPosition = (Mmatrix * vec4(position, 1.)).xyz;
            }`;
            const fsSource = `
            precision mediump float;
            varying vec3 vColor;
            varying vec2 vTexcoord;
            varying vec3 vNormal;
            varying vec3 v_worldPosition;
            uniform sampler2D sampler;
            uniform int u_useTexture;
            uniform vec3 u_lightPosition;

            void main(void) {
                vec3 normal = normalize(vNormal);
                vec3 lightDirection = normalize(u_lightPosition - v_worldPosition);
                float diff = max(dot(normal, lightDirection), 0.0);
                vec3 diffuse = diff * vec3(1.0, 1.0, 1.0);

                vec4 baseColor;
                if (u_useTexture == 1) {
                    baseColor = texture2D(sampler, vTexcoord);
                } else {
                    baseColor = vec4(vColor, 1.);
                }
                // Apply ambient light (0.9) + diffuse light (0.2)
                gl_FragColor = vec4(baseColor.rgb * (0.9 + diffuse * 0.2), baseColor.a);
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
                u_useTexture: gl.getUniformLocation(program, "u_useTexture"),
                u_lightPosition: gl.getUniformLocation(program, "u_lightPosition")
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

        /**
         * Creates a new shader program specifically for the skybox.
         * This shader is simple: it's not affected by lighting.
         */
        createSkyboxShaderProgram() {
            const gl = this.gl;
            const vsSource = `
            attribute vec3 position;
            varying vec3 v_texCoord;
            uniform mat4 Vmatrix;
            uniform mat4 Pmatrix;
            void main(void) {
                // Use the vertex position as the 3D texture coordinate
                v_texCoord = position;
                
                // Set the view matrix translation to 0,0,0
                // This makes the skybox "follow" the camera
                mat4 VmatrixNoTranslation = Vmatrix;
                VmatrixNoTranslation[3][0] = 0.0;
                VmatrixNoTranslation[3][1] = 0.0;
                VmatrixNoTranslation[3][2] = 0.0;
                
                gl_Position = Pmatrix * VmatrixNoTranslation * vec4(position, 1.0);
                
                // Force the skybox to be at the far clip plane (z/w = 1.0)
                // This ensures it's always drawn behind everything else.
                gl_Position = gl_Position.xyww;
            }`;

            const fsSource = `
            precision mediump float;
            varying vec3 v_texCoord;
            uniform samplerCube u_skybox;
            void main(void) {
                // Look up the color from the cubemap
                gl_FragColor = textureCube(u_skybox, normalize(v_texCoord));
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
                Pmatrix: gl.getUniformLocation(program, "Pmatrix"),
                Vmatrix: gl.getUniformLocation(program, "Vmatrix"),
                u_skybox: gl.getUniformLocation(program, "u_skybox")
            };

            gl.enableVertexAttribArray(locations.position);
            return { program, locations };
        }

        /**
         * Creates the vertex and face buffers for the skybox.
         */
        createSkyboxBuffers() {
            const gl = this.gl;
            const skyboxGeo = Geometry.generateSkyboxCube(1.0); // Size doesn't matter, we scale it

            const vertexBuffer = gl.createBuffer();
            gl.bindBuffer(gl.ARRAY_BUFFER, vertexBuffer);
            gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(skyboxGeo.vertices), gl.STATIC_DRAW);

            const facesBuffer = gl.createBuffer();
            gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, facesBuffer);
            gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(skyboxGeo.faces), gl.STATIC_DRAW);

            return {
                vertex: vertexBuffer,
                faces: facesBuffer,
                faces_length: skyboxGeo.faces.length
            };
        }

        /**
         * Loads 6 images into a single CUBE_MAP texture.
         */
        loadCubemapTexture(urls) {
            const gl = this.gl;
            const texture = gl.createTexture();
            gl.bindTexture(gl.TEXTURE_CUBE_MAP, texture);

            const targets = [
                gl.TEXTURE_CUBE_MAP_POSITIVE_X, gl.TEXTURE_CUBE_MAP_NEGATIVE_X,
                gl.TEXTURE_CUBE_MAP_POSITIVE_Y, gl.TEXTURE_CUBE_MAP_NEGATIVE_Y,
                gl.TEXTURE_CUBE_MAP_POSITIVE_Z, gl.TEXTURE_CUBE_MAP_NEGATIVE_Z
            ];

            let imagesLoaded = 0;

            for (let i = 0; i < 6; i++) {
                const image = new Image();
                image.crossOrigin = "anonymous"; // In case you load from another domain
                image.onload = () => {
                    gl.bindTexture(gl.TEXTURE_CUBE_MAP, texture);
                    gl.texImage2D(targets[i], 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, image);
                    imagesLoaded++;

                    if (imagesLoaded === 6) {
                        gl.generateMipmap(gl.TEXTURE_CUBE_MAP);
                        gl.texParameteri(gl.TEXTURE_CUBE_MAP, gl.TEXTURE_MIN_FILTER, gl.LINEAR_MIPMAP_LINEAR);
                    }
                };
                image.src = urls[i];
            }

            gl.texParameteri(gl.TEXTURE_CUBE_MAP, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
            gl.texParameteri(gl.TEXTURE_CUBE_MAP, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);

            return texture;
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

                // MODIFIED: Apply the rotation to the *entire environment*
                this.environment.modelMatrix = rotationMatrix;
            };
        }

        startRenderLoop() {
            const gl = this.gl;
            gl.enable(gl.DEPTH_TEST);
            gl.depthFunc(gl.LEQUAL);
            gl.clearColor(0.0, 0.0, 0.0, 0.0);
            gl.clearDepth(1.0);

            const render = (now) => {
                this.updateRotation();

                // --- CALCULATE ANIMATION VALUES ---
                const timeInSeconds = now * 0.001;
                const bodyBreathSpeed = 1.5;
                const bodyBreathAmount = 0.04;
                const bodyBreathOffset  = Math.sin(timeInSeconds * bodyBreathSpeed * Math.PI) * bodyBreathAmount;

                const eyeBreathSpeed = 1.5;
                const eyeBreathAmount = 0.02;
                const eyeBreathScale = 1.0 + Math.sin(timeInSeconds * eyeBreathSpeed * Math.PI) * eyeBreathAmount;

                const diskBreathSpeed = 1.0;
                const diskBreathAmount = 0.01;
                const diskBreathScale = 1.0 + Math.cos(timeInSeconds * diskBreathSpeed * Math.PI) * diskBreathAmount;

                const flapSpeed = 2.0;
                const flapAmount = 0.2; // Radians
                const flapAngle = Math.sin(timeInSeconds * flapSpeed * Math.PI) * flapAmount;

                const animationValues = {
                    body: bodyBreathOffset,
                    eye: eyeBreathScale,
                    disk: diskBreathScale,
                    flapAngle: flapAngle,
                };
                // --- END ANIMATION VALUES ---

                // --- UPDATE MODELS ---
                this.environment.updateAnimation(); // Updates ice island float

                // NEW: Update Piplup's animation
                // Note: Piplup's update function just takes 'time', not 'animationValues'
                this.piplup.updateAnimation(timeInSeconds * 7); // *1.5 to make it run a bit faster

                this.Prinplup.updateAnimation(animationValues); // Updates body/hand anims

                // Empoleon update animation
                this.empoleon.updateAnimation(animationValues);


                // --- RENDER ---
                gl.viewport(0, 0, this.canvas.width, this.canvas.height);
                gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

                // === 1. DRAW SKYBOX ===
                gl.depthMask(false); // Disable writing to the depth buffer
                gl.useProgram(this.skyboxShader.program);

                // NEW: Disable attributes from the main shader that we don't use
                gl.disableVertexAttribArray(this.shader.locations.color);
                gl.disableVertexAttribArray(this.shader.locations.texcoord);
                gl.disableVertexAttribArray(this.shader.locations.normal);
                // We leave 'position' enabled as both shaders use it

                // NEW: Create a view matrix for the skybox that includes the mouse rotation
                // this.environment.modelMatrix holds the rotation from the mouse
                const skyboxViewMatrix = LIBS.multiply(this.viewMatrix, this.environment.modelMatrix);

                // MODIFIED: Use the new skyboxViewMatrix instead of this.viewMatrix
                gl.uniformMatrix4fv(this.skyboxShader.locations.Pmatrix, false, this.projMatrix);
                gl.uniformMatrix4fv(this.skyboxShader.locations.Vmatrix, false, skyboxViewMatrix);

                gl.activeTexture(gl.TEXTURE0);
                gl.bindTexture(gl.TEXTURE_CUBE_MAP, this.skyboxTexture);
                gl.uniform1i(this.skyboxShader.locations.u_skybox, 0);

                // Bind buffers
                gl.bindBuffer(gl.ARRAY_BUFFER, this.skyboxBuffers.vertex);
                // Tell WebGL how to read the position-only buffer
                gl.vertexAttribPointer(this.skyboxShader.locations.position, 3, gl.FLOAT, false, 0, 0);

                gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, this.skyboxBuffers.faces);

                // Draw
                gl.drawElements(gl.TRIANGLES, this.skyboxBuffers.faces_length, gl.UNSIGNED_SHORT, 0);

                gl.depthMask(true); // Re-enable depth writing

                // === 2. DRAW MAIN SCENE ===
                gl.useProgram(this.shader.program); // Switch back to the main shader

                // NEW: Re-enable attributes for the main shader
                gl.enableVertexAttribArray(this.shader.locations.color);
                gl.enableVertexAttribArray(this.shader.locations.texcoord);
                gl.enableVertexAttribArray(this.shader.locations.normal);

                // Set main scene uniforms
                gl.uniformMatrix4fv(this.shader.locations.Pmatrix, false, this.projMatrix);
                gl.uniformMatrix4fv(this.shader.locations.Vmatrix, false, this.viewMatrix);
                gl.uniform3fv(this.shader.locations.u_lightPosition, [5, 15, 10]);

                // Get the final animated matrix of the ice island
                const iceParentMatrix = this.environment.getIceIslandWorldMatrix();

                // Draw the models
                this.Prinplup.draw(this.shader, iceParentMatrix);

                // Draw Piplup
                // First, combine the ice matrix with Piplup's static offset
                const piplupParentMatrix = LIBS.multiply(this.piplupModelMatrix, iceParentMatrix);
                // Then, draw Piplup using this final matrix
                this.piplup.draw(this.shader, piplupParentMatrix);

                // Draw Empoleon
                const empoleonParentMatrix = LIBS.multiply(this.empoleonModelMatrix, iceParentMatrix);
                this.empoleon.draw(this.shader, empoleonParentMatrix);

                this.environment.draw(this.shader);

                requestAnimationFrame(render);
            };
            render();
        }
    }

// --- START THE APPLICATION ---
    window.addEventListener('load', () => {
        // Make sure your canvas ID matches
        new Renderer('prinplup-canvas');
    });

})();