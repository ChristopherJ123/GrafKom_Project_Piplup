# Laporan Pengerjaan Proyek WebGL Grafika Komputer
<p>Kelompok 15: Piplup, Prinplup, Empoleon</p>

## Breakdown bentuk final Pokemon:
### Piplup (Chris)
Piplup dibuat menggunakan struktur hierarki `ModelNode`. Bagian-bagian utama seperti badan, kepala, mata, paruh, tangan (flipper), kaki, dan detail jubah/ekor dibuat menggunakan geometri dasar seperti `generateSphere`, `generateCircle`, dan `generateBeakHalf` (untuk paruh). Setiap bagian diposisikan relatif terhadap induknya (misalnya, mata adalah anak dari kepala, kepala adalah anak dari badan). Tekstur diterapkan pada kepala (piplup_head_texture.png). Struktur hierarkis ini memungkinkan animasi yang kompleks dan terkoordinasi.
### Prinplup (Michelle)
Prinplup awalnya dideskripsikan menggunakan 3 elipsoid tumpuk untuk badan. Namun, implementasi dalam kode prinplup2_lighting.js menggunakan teknik `generateLathe` untuk membuat bentuk object terkait dari control points body (body_profile), menghasilkan bentuk yang menyerupai profil tubuh Prinplup yang menyerupai vas bunga. Bagian lain seperti kepala (termasuk disk kuning), mata (`generateAngryEye`), paruh (`generateBeak` dan `generateCone`), tangan (`generateSphere`), kaki (`generateSphere`), dan sirip belakang (`generateTaperedShapeFromSpline`) ditambahkan sebagai node terpisah. Tekstur diterapkan pada badan (prinplup_texture_bare.png) dan tangan (prinplup_hand_texture.png).
### Empoleon (Mario)
Empoleon juga menggunakan `generateLathe` untuk object terkait berdasarkan control points body (body_profile). Bentuk-bentuk kompleks seperti trisula di kepala menggunakan kombinasi dari (`generateTriangularPrism`) untuk bagian tiang mahkotanya dan 2 generateCone (sisi atas dan bawah) untuk membentuk diamond di ujung mahkotanya. Pada paruhnya menggunakan (`generateTriBeak`) yang bisa dicustom panjang lebar dan tingginya. Kemudian detail jubah/tangan menggunakan kombinasi (`generateTaperedShapeFromSpline`) untuk jubahnya, (`generateHalfEllipsoid`) untuk bagian outer hand, (`generateCone`) untuk bagian tajam lengan dan juga jari, (`generateSphere`) untuk lengannya. Kaki menggunakan (`generateHalfHyperboloid untuk`) bagian 'Legs' dan (`generateSphere`) untuk bagian 'Feet'. 2 Bulatan di punggung berwarna kuning menggunakan (`generateSphere`), kemudian untuk ekor menggunakan (`generateTaperedShapeFromSpline`). Tekstur diterapkan pada object terkait (empoleon_texture.png).

<div style="margin-top: 5rem;">

## Breakdown environment: 
Elemen environment yang digunakan adalah pulau es (iceberg) tempat object Pokemon berdiri dan bidang air di sekelilingnya.
- Pulau Es:
    <p>Dibuat menggunakan fungsi custom <code>generateIrregularExtrudedPolygon</code>. Fungsi ini menghasilkan mesh 3D berbentuk poligon (segi-N, contohnya 6 untuk hexagon) yang diekstrusi (diberi ketebalan). Posisi titik-titik sudut poligon diacak sedikit menggunakan angka random dan <code>irregularityFactor</code> untuk membuatnya tidak simetris agar menyerupai bongkahan es alami.</p>
- Water plane:
    <p>Dibuat menggunakan fungsi <code>generateWaterPlane</code>. Fungsi ini membuat bidang (plane) namun posisi Y (ketinggian) setiap vertex diubah menggunakan kombinasi fungsi sin dan cos berdasarkan koordinat X dan Z yang menghasilkan efek gelombang statis. Texture diterapkan pada object terkait (water-texture.png).</p>
    <i style="margin-bottom: 2rem;">Kalau mau diganti hanya warna solid, tinggal hapus <code>waterTexture</code> dari <code>waterNode</code></i>
- Background (Sky Box): 
    <p>Dibuat dengan SkyBox yang memanfaatkan keenam faces dari cube dengan image masing-masing</p>

## Penjelasan bentuk selain quadric: 
- `generateLathe`: Digunakan untuk membuat object terkait Prinplup dan Empoleon dengan memutar profil 2D di sekitar sumbu Y.

- `generateCone`: Digunakan untuk bagian atas paruh Prinplup dan detail trisula pada tangan Empoleon.

- `generateBeak` / `generateBeakHalf`: Fungsi custom untuk membuat bentuk paruh yang meruncing dan sedikit melengkung (seperti paraboloid eliptik).

- `generateAngryEye`: Fungsi custom (variasi dari sphere/ellipsoid) untuk membuat bentuk mata Prinplup yang menyipit dengan membatasi sapuan vertikal saat generasi vertex.

- `generateCircle`: Digunakan untuk membuat detail lingkaran putih datar di dada Piplup dan Prinplup.

- `generateTubeFromSpline`: Digunakan untuk membuat detail jubah/ekor Piplup yang melengkung mengikuti jalur spline.

- `generateTaperedShapeFromSpline`: Digunakan untuk membuat sirip ekor Prinplup dan Empoleon yang bentuknya meruncing dari pangkal ke ujung mengikuti jalur spline.

- `generateIrregularExtrudedPolygon`: Fungsi custom untuk pulau es, membuat mesh 3D dari poligon 2D yang tidak beraturan.

- `generateWaterPlane`: Fungsi custom untuk bidang air, menghasilkan mesh plane dengan permukaan bergelombang.

- `generateHalfEllipsoid`: Digunakan untuk detail tangan Empoleon.

- `generateHalfHyperboloid`: Digunakan untuk kaki Empoleon.

## Breakdown animasi: 
### Piplup (Struktur Hierarkis ModelNode):
- Lari (Bobbing & Squash/Stretch): Badan utama (`bodyNode` dan `bodyGeometry`) digerakkan naik-turun (`translateY`) dan diskalakan (`scaleY`, `scaleX`, `scaleZ`) untuk efek memantul dan melentur.
- Ayunan Tangan (Flipper Flapping): Tangan (`leftFlipper`, `rightFlipper`) dirotasi pada sumbu Z  dengan titik pivot di pangkal lengan menggunakan komposisi translasi-rotasi-translasi.
- Ayunan Kaki: Kaki (`leftLeg`, `rightLeg`) dirotasi pada sumbu X (`rotateX`) secara bergantian.
- Paruh Berbicara: Bagian atas dan bawah paruh (`topBeak`, `bottomBeak`) dirotasi pada sumbu X (`rotateX`) dengan titik pivot di pangkal paruh.
- Animasi di-update di `Piplup.updateAnimation(time)` dan diterapkan pada `localMatrix` setiap `ModelNode` terkait. Matrix akhir dihitung secara rekursif melalui `updateWorldMatrix`.
### Prinplup (Pola animValues & PrinplupPart):
- Bernapas (Badan): Nilai `bodyBreathValue` (hasil dari sin) dihitung di `Prinplup.draw()` dan diteruskan ke `PrinplupPart.draw()`. Bagian dengan `animationType: 'bodyBreathe'` akan menerapkan translasi Y menggunakan nilai ini.
- Bernapas (Mata): Nilai `eyeBreathValue` (hasil dari 1.0 + sin) dihitung dan diteruskan. Bagian mata (`animationType: 'eyeBreathe'`) menerapkan Y scale.
- Bernapas (Disk): Nilai `diskBreathValue` (hasil dari 1.0 + cos) dihitung dan diteruskan. Bagian disk (`animationType: 'diskBreathe'`) menerapkan uniform scale.
- Tangan Clapping: Nilai `clapRotationValue` (hasil dari sin absolute) dihitung dan diteruskan. Bagian tangan kiri (`animationType: 'handClapLeft'`) menerapkan rotasi Z positif, tangan kanan (`animationType: 'handClapRight'`) menerapkan rotasi Z negatif.
- Animasi diterapkan di `PrinplupPart.draw()` dengan mengalikan matrix posisi asli (`this.modelMatrix`) dengan matrix animasi lokal (`animMatrix`) yang dibuat berdasarkan `animValues`.
### Environment (Pulau Es Mengambang):
- Mengambang: Nilai floatY (hasil dari sin) dihitung di `Environment.draw()`. Matrix translasi Y (`iceAnimMatrix`) dibuat dan disimpan. Matrix ini juga digunakan oleh `PrinplupPart.draw()` untuk bagian pulau es (`animationType: 'floatingIce'`).
- Parenting: Renderer mengambil `iceAnimMatrix` dari Environment (melalui `getIceIslandAnimMatrix()`) dan meneruskannya sebagai `parentMatrix` saat memanggil `Prinplup.draw()`, sehingga Prinplup ikut bergerak naik-turun bersama pulau es.

</div>


<div style="margin-top: 9rem;">

## Kriteria Penilaian: 
Pastikan semua kriteria penilaian telah dijelaskan dan ditunjukkan dalam laporan ini.

## Penjelasan tambahan (jika ada): 
Jika kalian menggunakan fitur tambahan seperti lighting, opacity, tekstur, efek kamera, dan sebagainya, jelaskan juga secara singkat penggunaannya.
</div>
