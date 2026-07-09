#!/bin/bash

# clean up previous output files
# rm -f Tracking/test/testTrackFinder/out_sim_edm4hep.root
# rm -f Tracking/test/testTrackFinder/out_tracks.root

MODEL_PATH=$1

XML_FILE=$K4GEO/FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml
STEERING_FILE=SteeringFile_IDEA_o1_v03.py
TBETA=0.6
TD=0.3

# curl -o $STEERING_FILE https://raw.githubusercontent.com/key4hep/k4geo/master/example/SteeringFile_IDEA_o1_v03.py

# ddsim --steeringFile $STEERING_FILE \
#       --compactFile  $XML_FILE \
#       -G --gun.distribution uniform --gun.particle mu- \
#       --random.seed 42 \
#       --numberOfEvents 1 \
#       --outputFile out_sim_edm4hep.root

k4run runTestTrackFinder.py --inputFile out_sim_edm4hep.root --outputFile out_tracks.root --modelPath $MODEL_PATH --tbeta $TBETA --td $TD


# GGTFTrackFinder       INFO Hit 0 input features: (13.6085, 3.17897, -24.5444, 1, 0, 0, 0)
# GGTFTrackFinder       INFO Hit 1 input features: (23.725, 5.55789, -42.7971, 1, 0, 0, 0)
# GGTFTrackFinder       INFO Hit 2 input features: (34.7349, 8.16486, -62.6681, 1, 0, 0, 0)
# GGTFTrackFinder       INFO Hit 3 input features: (170.558, 41.6827, -308.415, 1, 0, 0, 0)
# GGTFTrackFinder       INFO Hit 4 input features: (1253.69, 399.473, -2314.17, 1, 0, 0, 0)
# GGTFTrackFinder       INFO Hit 5 input features: (1264.24, 404.16, -2334.17, 1, 0, 0, 0)
# GGTFTrackFinder       INFO Hit 6 input features: (1073.49, 329.163, -1965.41, 0, 0, -1.80788, -0.262728)
# GGTFTrackFinder       INFO Hit 7 input features: (1015.1, 319.248, -1890.84, 0, 0, -9.26804, -1.27739)
# GGTFTrackFinder       INFO Hit 8 input features: (987.174, 303.252, -1785.82, 0, 0, -10.6432, -1.42753)
# GGTFTrackFinder       INFO Hit 9 input features: (1057.88, 328.136, -1938.82, 0, 0, -4.22954, 0.51662)
# GGTFTrackFinder       INFO Hit 10 input features: (971.74, 297.34, -1765.52, 0, 0, -4.12469, 0.471996)
# GGTFTrackFinder       INFO Hit 11 input features: (959.125, 285.267, -1753.08, 0, 0, -0.465818, -0.060687)
# GGTFTrackFinder       INFO Hit 12 input features: (931.518, 283.203, -1767.36, 0, 0, -8.93214, -1.13032)
# GGTFTrackFinder       INFO Hit 13 input features: (877.866, 248.94, -1546.79, 0, 0, -12.1899, -1.45343)
# GGTFTrackFinder       INFO Hit 14 input features: (860.023, 257.17, -1575.13, 0, 0, -6.61027, 0.683773)
# GGTFTrackFinder       INFO Hit 15 input features: (845.46, 254.886, -1566.45, 0, 0, -7.51848, -0.866534)
# GGTFTrackFinder       INFO Hit 16 input features: (1001.2, 313.125, -1843.84, 0, 0, -11.5131, 1.3452)
# GGTFTrackFinder       INFO Hit 17 input features: (915.614, 278.831, -1681.58, 0, 0, -8.88032, 0.966775)
# GGTFTrackFinder       INFO Hit 18 input features: (834.226, 239.353, -1536.64, 0, 0, -6.85242, 0.691967)
# GGTFTrackFinder       INFO Hit 19 input features: (818.79, 236.644, -1506.64, 0, 0, -0.636791, -0.0710297)
# GGTFTrackFinder       INFO Hit 20 input features: (804.053, 238.316, -1471.28, 0, 0, -7.19416, 0.703086)
# GGTFTrackFinder       INFO Hit 21 input features: (764.788, 217.016, -1363.56, 0, 0, -3.76099, -0.392077)
# GGTFTrackFinder       INFO Hit 22 input features: (710.524, 196.476, -1268.77, 0, 0, -11.9557, -1.15654)
# GGTFTrackFinder       INFO Hit 23 input features: (697.727, 191.664, -1299.38, 0, 0, -8.4237, 0.728544)
# GGTFTrackFinder       INFO Hit 24 input features: (680.807, 201.01, -1276.83, 0, 0, -8.01698, -0.745611)
# GGTFTrackFinder       INFO Hit 25 input features: (667.795, 191.126, -1187.51, 0, 0, -2.51229, 0.209305)
# GGTFTrackFinder       INFO Hit 26 input features: (656.093, 182.467, -1205.46, 0, 0, -2.88126, -0.257823)
# GGTFTrackFinder       INFO Hit 27 input features: (642.907, 175.448, -1187.43, 0, 0, -5.76097, 0.463305)
# GGTFTrackFinder       INFO Hit 28 input features: (612.949, 174.035, -1114.5, 0, 0, -6.43901, 0.495724)
# GGTFTrackFinder       INFO Hit 29 input features: (587.789, 162.21, -1037.45, 0, 0, -11.5611, 0.860153)
# GGTFTrackFinder       INFO Hit 30 input features: (561.592, 153.32, -1000.36, 0, 0, -1.90803, 0.135898)
# GGTFTrackFinder       INFO Hit 31 input features: (548.28, 149.738, -1013.45, 0, 0, -0.882746, -0.0660228)
# GGTFTrackFinder       INFO Hit 32 input features: (534.595, 148.119, -1022.72, 0, 0, -8.24696, 0.558394)
# GGTFTrackFinder       INFO Hit 33 input features: (521.045, 138.361, -934.601, 0, 0, -4.53609, -0.321958)
# GGTFTrackFinder       INFO Hit 34 input features: (598.896, 174.557, -1075.49, 0, 0, -12.3607, -1.01097)
# GGTFTrackFinder       INFO Hit 35 input features: (508.764, 134.198, -945.648, 0, 0, -11.7267, 0.76156)
# GGTFTrackFinder       INFO Hit 36 input features: (1033.15, 309.763, -1918.58, 0, 0, -6.4787, 0.777061)
# GGTFTrackFinder       INFO Hit 37 input features: (947.942, 276.724, -1765.55, 0, 0, -7.30785, 0.819068)
# GGTFTrackFinder       INFO Hit 38 input features: (494.511, 133.541, -929.035, 0, 0, -9.51015, -0.640566)
# GGTFTrackFinder       INFO Hit 39 input features: (480.79, 130.97, -861.842, 0, 0, -8.05285, 0.496292)
# GGTFTrackFinder       INFO Hit 40 input features: (466.478, 137.134, -874.013, 0, 0, -12.349, -0.786353)
# GGTFTrackFinder       INFO Hit 41 input features: (456.684, 121.771, -880.126, 0, 0, -4.75701, 0.278215)
# GGTFTrackFinder       INFO Hit 42 input features: (889.099, 262.475, -1611.45, 0, 0, -2.94689, 0.314139)
# GGTFTrackFinder       INFO Hit 43 input features: (440.932, 130.273, -823.974, 0, 0, -13.0403, -0.784647)
# GGTFTrackFinder       INFO Hit 44 input features: (1045.8, 322.821, -1956.68, 0, 0, -8.43163, -1.19365)
# GGTFTrackFinder       INFO Hit 45 input features: (427.343, 120.368, -768.445, 0, 0, -8.24359, 0.453949)
# GGTFTrackFinder       INFO Hit 46 input features: (781.4, 218.499, -1463.28, 0, 0, -6.67765, 0.636874)
# GGTFTrackFinder       INFO Hit 47 input features: (413.918, 114.385, -724.874, 0, 0, -8.13438, -0.458441)
# GGTFTrackFinder       INFO Hit 48 input features: (792.896, 227.26, -1453.34, 0, 0, -8.50021, -0.917615)
# GGTFTrackFinder       INFO Hit 49 input features: (403.236, 103.285, -747.757, 0, 0, -12.9455, 0.675769)
# GGTFTrackFinder       INFO Hit 50 input features: (901.75, 275.372, -1658.62, 0, 0, -9.69052, -1.19)
# GGTFTrackFinder       INFO Hit 51 input features: (387.8, 104.66, -712.6, 0, 0, -3.97226, -0.209686)
# GGTFTrackFinder       INFO Hit 52 input features: (376.437, 97.2791, -674.189, 0, 0, -6.02651, 0.294659)
# GGTFTrackFinder       INFO Hit 53 input features: (363.788, 92.9042, -644.242, 0, 0, -0.241192, -0.0119253)
# GGTFTrackFinder       INFO Hit 54 input features: (351.819, 90.7349, -626.002, 0, 0, -0.482232, 0.0220842)
# GGTFTrackFinder       INFO Hit 0 model output: (-1.7047, 1.02814, 1.11461, -3.50257) cluster ID: 1
# GGTFTrackFinder       INFO Hit 1 model output: (-1.71787, 1.04335, 1.12203, -4.00938) cluster ID: 1
# GGTFTrackFinder       INFO Hit 2 model output: (-1.71728, 1.04, 1.12245, -3.62614) cluster ID: 1
# GGTFTrackFinder       INFO Hit 3 model output: (-1.71823, 1.04758, 1.10757, -6.4146) cluster ID: 1
# GGTFTrackFinder       INFO Hit 4 model output: (-1.75054, 1.05659, 1.1312, -3.13408) cluster ID: 1
# GGTFTrackFinder       INFO Hit 5 model output: (-1.73067, 1.04918, 1.12212, -2.87726) cluster ID: 1
# GGTFTrackFinder       INFO Hit 6 model output: (-1.76508, 1.03758, 1.13145, 3.69901) cluster ID: 1
# GGTFTrackFinder       INFO Hit 7 model output: (-1.76427, 1.06473, 1.13454, -3.99672) cluster ID: 1
# GGTFTrackFinder       INFO Hit 8 model output: (-1.75156, 1.06965, 1.13015, -5.21521) cluster ID: 1
# GGTFTrackFinder       INFO Hit 9 model output: (-1.76873, 1.06873, 1.13309, -2.79991) cluster ID: 1
# GGTFTrackFinder       INFO Hit 10 model output: (-1.75823, 1.06503, 1.12986, -4.55688) cluster ID: 1
# GGTFTrackFinder       INFO Hit 11 model output: (-1.76079, 1.06059, 1.12993, -4.30028) cluster ID: 1
# GGTFTrackFinder       INFO Hit 12 model output: (-1.74408, 1.0481, 1.12165, -3.87641) cluster ID: 1
# GGTFTrackFinder       INFO Hit 13 model output: (-1.73792, 1.06158, 1.11728, -4.27492) cluster ID: 1
# GGTFTrackFinder       INFO Hit 14 model output: (-1.74677, 1.0645, 1.12715, -4.86508) cluster ID: 1
# GGTFTrackFinder       INFO Hit 15 model output: (-1.7385, 1.05509, 1.12389, -4.68134) cluster ID: 1
# GGTFTrackFinder       INFO Hit 16 model output: (-1.75428, 1.06553, 1.13308, -4.72253) cluster ID: 1
# GGTFTrackFinder       INFO Hit 17 model output: (-1.75296, 1.06611, 1.12993, -4.5253) cluster ID: 1
# GGTFTrackFinder       INFO Hit 18 model output: (-1.73508, 1.04404, 1.11855, -5.02076) cluster ID: 1
# GGTFTrackFinder       INFO Hit 19 model output: (-1.74599, 1.05795, 1.12318, -4.5693) cluster ID: 1
# GGTFTrackFinder       INFO Hit 20 model output: (-1.73547, 1.05, 1.12305, -4.91888) cluster ID: 1
# GGTFTrackFinder       INFO Hit 21 model output: (-1.73854, 1.06389, 1.11852, -4.95696) cluster ID: 1
# GGTFTrackFinder       INFO Hit 22 model output: (-1.73255, 1.04818, 1.11494, -5.49565) cluster ID: 1
# GGTFTrackFinder       INFO Hit 23 model output: (-1.73432, 1.04433, 1.11748, -4.3228) cluster ID: 1
# GGTFTrackFinder       INFO Hit 24 model output: (-1.7274, 1.05163, 1.11717, -4.95013) cluster ID: 1
# GGTFTrackFinder       INFO Hit 25 model output: (-1.73395, 1.06154, 1.12, -5.22381) cluster ID: 1
# GGTFTrackFinder       INFO Hit 26 model output: (-1.73379, 1.05477, 1.11953, -4.89871) cluster ID: 1
# GGTFTrackFinder       INFO Hit 27 model output: (-1.7206, 1.03974, 1.11288, -5.30822) cluster ID: 1
# GGTFTrackFinder       INFO Hit 28 model output: (-1.7296, 1.05639, 1.12143, -4.69312) cluster ID: 1
# GGTFTrackFinder       INFO Hit 29 model output: (-1.72366, 1.05351, 1.11559, -5.38911) cluster ID: 1
# GGTFTrackFinder       INFO Hit 30 model output: (-1.72837, 1.05392, 1.11616, -5.15473) cluster ID: 1
# GGTFTrackFinder       INFO Hit 31 model output: (-1.72659, 1.04512, 1.11306, -4.54574) cluster ID: 1
# GGTFTrackFinder       INFO Hit 32 model output: (-1.72182, 1.03263, 1.11003, -5.82318) cluster ID: 1
# GGTFTrackFinder       INFO Hit 33 model output: (-1.71996, 1.04078, 1.1112, -5.13186) cluster ID: 1
# GGTFTrackFinder       INFO Hit 34 model output: (-1.7278, 1.06005, 1.11934, -5.53942) cluster ID: 1
# GGTFTrackFinder       INFO Hit 35 model output: (-1.72406, 1.03687, 1.11045, -4.91466) cluster ID: 1
# GGTFTrackFinder       INFO Hit 36 model output: (-1.765, 1.06134, 1.13335, -3.70906) cluster ID: 1
# GGTFTrackFinder       INFO Hit 37 model output: (-1.7471, 1.04685, 1.12247, -3.94589) cluster ID: 1
# GGTFTrackFinder       INFO Hit 38 model output: (-1.71645, 1.03149, 1.10915, -5.03734) cluster ID: 1
# GGTFTrackFinder       INFO Hit 39 model output: (-1.71828, 1.04794, 1.11301, -5.40082) cluster ID: 1
# GGTFTrackFinder       INFO Hit 40 model output: (-1.72835, 1.05133, 1.1167, -5.65002) cluster ID: 1
# GGTFTrackFinder       INFO Hit 41 model output: (-1.72454, 1.05316, 1.11225, -5.12349) cluster ID: 1
# GGTFTrackFinder       INFO Hit 42 model output: (-1.75375, 1.07126, 1.12777, -4.96368) cluster ID: 1
# GGTFTrackFinder       INFO Hit 43 model output: (-1.73366, 1.05407, 1.11503, -5.37915) cluster ID: 1
# GGTFTrackFinder       INFO Hit 44 model output: (-1.75839, 1.0654, 1.13273, -4.28482) cluster ID: 1
# GGTFTrackFinder       INFO Hit 45 model output: (-1.72784, 1.05288, 1.11502, -5.24741) cluster ID: 1
# GGTFTrackFinder       INFO Hit 46 model output: (-1.73957, 1.0404, 1.11791, -4.41167) cluster ID: 1
# GGTFTrackFinder       INFO Hit 47 model output: (-1.72491, 1.04963, 1.11702, -5.73747) cluster ID: 1
# GGTFTrackFinder       INFO Hit 48 model output: (-1.73285, 1.04666, 1.11692, -5.04824) cluster ID: 1
# GGTFTrackFinder       INFO Hit 49 model output: (-1.73456, 1.04889, 1.1155, -5.31015) cluster ID: 1
# GGTFTrackFinder       INFO Hit 50 model output: (-1.74189, 1.0613, 1.12694, -4.78299) cluster ID: 1
# GGTFTrackFinder       INFO Hit 51 model output: (-1.73092, 1.05012, 1.11324, -4.99032) cluster ID: 1
# GGTFTrackFinder       INFO Hit 52 model output: (-1.73363, 1.05315, 1.11754, -6.16194) cluster ID: 1
# GGTFTrackFinder       INFO Hit 53 model output: (-1.73581, 1.05536, 1.1116, -5.22282) cluster ID: 1
# GGTFTrackFinder       INFO Hit 54 model output: (-1.74729, 1.05894, 1.12063, -5.69039) cluster ID: 1