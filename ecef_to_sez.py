# ecef_to_sez.py
#
# Usage: python3 ecef_to_sez.py o_x_km o_y_km o_z_km x_km y_km z_km
#  Converts sez to ECEF vector components
# Parameters:
#  o_x_deg: x ECEF vector origin of the SEZ frame
#  o_y_deg: y ECEF vector origin of the SEZ frame
#  o_z_km: z ECEF vector origin of the SEZ frame
#  x_km: x vector of ECEF poisition
#  y_km: y vector of ECEF poisition
#  z_km: z vector of ECEF poisition
# Output:
#  Prints the SEZ fram in km s (km), e (km), and z (km)
#
# Written by Yonghwa Kim
# Other contributors: None

# import Python modules
import math # math module
import sys  # argv

# "constants"
R_E_KM = 6378.1363
E_E    = 0.081819221456

# helper functions

## calculated denominator
def calc_denom(ecc, lat_rad):
  return math.sqrt(1.0-(ecc**2)*(math.sin(lat_rad)**2))

# initialize script arguments
o_x_km = float('nan')
o_y_km = float('nan') 
o_z_km = float('nan') 
x_km = float('nan') 
y_km = float('nan') 
z_km = float('nan') 

# parse script arguments
if len(sys.argv)==7:
  o_x_km = float(sys.argv[1])
  o_y_km = float(sys.argv[2])
  o_z_km = float(sys.argv[3])
  x_km = float(sys.argv[4])
  y_km = float(sys.argv[5])
  z_km = float(sys.argv[6])
else:
  print(\
   'Usage: '\
   'python3 ecef_to_sez.py o_x_km o_y_km o_z_km x_km y_km z_km'\
  )
  exit()

# Determine the ECEF vector
r_x_km = o_x_km
r_y_km = o_y_km
r_z_km = o_z_km

#ECEF = [o_x_km-x_km,o_y_km-y_km,o_z_km-z_km]
ECEF = [x_km-o_x_km,y_km-o_y_km,z_km-o_z_km]

# calculate longitude
lon_rad = math.atan2(r_y_km,r_x_km)
lon_deg = lon_rad*180.0/math.pi

# initialize lat_rad, r_lon_km, r_z_km
lat_rad = math.asin(r_z_km/math.sqrt(r_x_km**2+r_y_km**2+r_z_km**2))
r_lon_km = math.sqrt(r_x_km**2+r_y_km**2)
prev_lat_rad = float('nan')

# iteratively find latitude
c_E = float('nan')
count = 0
while (math.isnan(prev_lat_rad) or abs(lat_rad-prev_lat_rad)>10e-7) and count<5:
  denom = calc_denom(E_E,lat_rad)
  c_E = R_E_KM/denom
  prev_lat_rad = lat_rad
  lat_rad = math.atan((r_z_km+c_E*(E_E**2)*math.sin(lat_rad))/r_lon_km)
  count = count+1
  
# calculate hae
hae_km = r_lon_km/math.cos(lat_rad)-c_E

# Inverse rotations calculations
Ry = [[math.sin(lat_rad), 0, -math.cos(lat_rad)], [0, 1, 0],[math.cos(lat_rad), 0, math.sin(lat_rad)]]
Rz = [[math.cos(lon_rad), math.sin(lon_rad), 0], [-math.sin(lon_rad), math.cos(lon_rad), 0], [0, 0, 1]]

# Function to perform matrix multiplication
def matrix_multiply(A, B):
    # A is m x n, B is n x p
    result = [[0 for _ in range(len(B[0]))] for _ in range(len(A))]
    for i in range(len(A)):
        for j in range(len(B[0])):
            for k in range(len(B)):
                result[i][j] += A[i][k] * B[k][j]
    return result

# Calculate Inverse_Rot
Inverse_Rot = matrix_multiply(Ry, Rz)

# Function to multiply a matrix by a vector
def matrix_vector_multiply(matrix, vector):
    result = [0 for _ in range(len(matrix))]
    for i in range(len(matrix)):
        for j in range(len(vector)):
            result[i] += matrix[i][j] * vector[j]
    return result

# Calculate SEZ
SEZ = matrix_vector_multiply(Inverse_Rot, ECEF)

s_km = SEZ[0]
e_km = SEZ[1]
z_km = SEZ[2]

#print(SEZ)
print(s_km)
print(e_km)
print(z_km)