#!/usr/bin/python

import sys
import os
import subprocess
from scipy.special import lambertw
from math import *

# from math import sqrt,ceil,pi


split_tol = 0.5
# split_tol = 1.0/sqrt(2.0)


def search_val(text, name, pre_delim, post_delim):
  if text == "":
    return 0.0
  else:
    matched_lines = [line for line in text if name in line]
    if post_delim == "":
      return float(str(matched_lines[0]).split(pre_delim)[1])
    else:
      return float(str(matched_lines[0]).split(pre_delim)[1].split(post_delim)[0])

def even(N):
  return int(2 * ceil(N / 2.0))

def oversampled_gridsize(M, m, flag):
  # if flag is a float it is interpreted as oversampling factor
  if flag != int(flag):
    return even(M*flag)

  # if flag is an integer it determines the way of constant oversampling
  if flag == -1:
    return M + 2*m+2
  if flag == -2:
    return 2*M
  else:
    return M + flag

# overwrite default parameters, if command line arguments are given
def set_str(nr, default):
  if len(sys.argv) > nr:
    if sys.argv[nr] >= 0:
      return str(sys.argv[nr])
  return default

def set_int(nr, default):
  if len(sys.argv) > nr:
    if sys.argv[nr] >= 0:
      return int(sys.argv[nr])
  return default

def set_float(nr, default):
  if len(sys.argv) > nr:
    if sys.argv[nr] > 0:
      return float(sys.argv[nr])
  return default

def set_number(nr, default):
  if len(sys.argv) > nr:
    if sys.argv[nr] > 0:
      s = sys.argv[nr]
      return float(s) if '.' in s else int(s)
  return default

def inv_kolper_potential_alpha(rc, L, Q2, tol):
  lamb = lambertw( 1.0 / tol * sqrt(Q2 * rc / L**3) )
  return 1.0/rc * sqrt(lamb.real)
def inv_kolper_field_alpha(rc, L, Q2, tol):
  l = log( 2.0 / tol * sqrt( Q2 / rc / L**3) )
  return 1.0/rc * sqrt(l)

def inv_kolper_potential_kc(alpha, L, Q2, tol):
  lamb = lambertw( 4.0 / 3.0 / L**2 * pow( sqrt(2.0 * Q2 / pi / alpha) / tol, 4.0/3.0) )
  return sqrt(3.0) / 2.0 * alpha * L / pi * sqrt(lamb.real)
def inv_kolper_field_kc(alpha, L, Q2, tol):
  lamb = lambertw(pow( sqrt(Q2 * alpha) / sqrt(pi * L**3) * 8.0 / tol , 4.0))
  return 0.5 * alpha * L / pi * sqrt(lamb.real)

def inv_kolper_alpha(rc, L, Q2, N, tol, tol_type):
  if tol_type == 1:
    return inv_kolper_potential_alpha(rc, L, Q2, tol)
  elif tol_type == 2:
    return inv_kolper_potential_alpha(rc, L, Q2, tol * sqrt(N/Q2))
  elif tol_type == 3:
    return inv_kolper_field_alpha(rc, L, Q2, tol)
  elif tol_type == 4:
    return inv_kolper_field_alpha(rc, L, Q2, tol * sqrt(N/Q2))
  else:
    print("ERROR: tol-type out of range")
    sys.exit()

def inv_kolper_kc(alpha, L, Q2, N, tol, tol_type):
  if tol_type == 1:
    return inv_kolper_potential_kc(alpha, L, Q2, tol)
  elif tol_type == 2:
    return inv_kolper_potential_kc(alpha, L, Q2, tol * sqrt(N/Q2))
  elif tol_type == 3:
    return inv_kolper_field_kc(alpha, L, Q2, tol)
  elif tol_type == 4:
    return inv_kolper_field_kc(alpha, L, Q2, tol * sqrt(N/Q2))
  else:
    print("ERROR: tol-type out of range")
    sys.exit()

def kolper_near_potential(rc, L, Q2, alpha):
  arg = pow(alpha*rc, 2.0)
  return sqrt(Q2 * rc / L**3) * exp(-arg) / arg
def kolper_near_field(rc, L, Q2, alpha):
  arg = pow(alpha*rc, 2.0)
  return 2.0 * sqrt(Q2) / sqrt(rc * pow(L,3.0)) * exp(-arg)

def kolper_far_potential(kc, L, Q2, alpha):
  arg = pow(pi * kc / (alpha * L), 2.0)
  return 2.0 * sqrt(Q2) * sqrt(0.5 / pow(kc,3.0)) * alpha/(pi*pi) * exp(-arg)
#   return 2.0 * sqrt(Q2/N * (alpha*L/kc)**2.0 + Q2 / 2.0 / pow(kc,3.0)) * alpha/(pi*pi) * exp(-arg)
def kolper_far_field(kc, L, Q2, alpha):
  arg = pow(pi * kc / (alpha * L), 2.0)
  # estimate by Franzisk Nestler
#   return 2.0 * alpha / L * sqrt(Q2) / sqrt(pi * kc) * exp(-arg)
  # estimate by Kolafa/Perram
  return 2.0 * alpha / (pi * L) * sqrt(8 * Q2 / kc) * exp(-arg)

def kolper_near_energy(rc, L, Q2, alpha, N):
  return sqrt(Q2/N) * kolper_near_potential(rc, L, Q2, alpha)
def kolper_far_energy(kc, L, Q2, alpha, N):
  return sqrt(Q2/N) * kolper_far_potential(kc, L, Q2, alpha)
def kolper_near_force(rc, L, Q2, alpha, N):
  return sqrt(Q2/N) * kolper_near_field(rc, L, Q2, alpha)
def kolper_far_force(kc, L, Q2, alpha, N):
  return sqrt(Q2/N) * kolper_far_field(kc, L, Q2, alpha)

def kolper_near_total_energy(rc, L, Q2, alpha):
#   if rc < L:
#     const = sqrt(pi) / L**4.0 / alpha**3.0 * Q2 * exp(-(alpha*L)**2.0)
#   else:
#     const = sqrt(pi) / rc / L**3.0 / alpha**3.0 * Q2 * exp(-(alpha*rc)**2.0)
#   cor = Q2 * energy_correction_near(rc, L, alpha)
#   print("Kolafa-Perram NEAR correction = "+ str('%.2e' % const))
#   print("Numerical NEAR correction     = "+ str('%.2e' % cor))
  return sqrt(Q2/2.0) * kolper_near_potential(rc, L, Q2, alpha) / 2.0
def kolper_far_total_energy(kc, L, Q2, alpha):
#   arg = pow(pi * kc / (alpha * L), 2.0)
#   const = Q2 * alpha**2.0 * L / pi**2.0 / kc * exp(-arg)
#   cor   = Q2 * energy_correction_far(kc, L, alpha)
#   print("Kolafa-Perram FAR correction = "+ str('%.2e' % const))
#   print("Numerical FAR correction     = "+ str('%.2e' % cor))
#   
#   rf20 = energy_correction_far(kc, L, alpha) 
#   arg = pow(pi * kc / (alpha * L), 2.0)
#   rEf2 = alpha/(pi*pi) * pow(kc, -1.5) * exp(-arg)
# 
#   f20 = rf20**2
#   Ef2 = rEf2**2
# 
#   Ef2 = alpha**2 / (2.0 * kc**3 * pi**4) * exp(-2.0*arg)
# 
#   h1 = (Q2*Q2-Q4) * f20
#   h2 =  Q2*Q2 * Ef2 
# 
#   print("h1 = "+ str("%.2e" % h1) + ", h2 = "+ str("%.e" % h2))
#   U = sqrt( h1 + h2 )

  Unc  = sqrt(Q2/2.0) * kolper_far_potential(kc, L, Q2, alpha) / 2.0
#   Uold = cor + Unc
# 
#   print("total energy non-corrected  : "+ str("%.2e" % Unc) )
#   print("total energy corrected (old): "+ str("%.2e" % Uold) )
#   print("total energy corrected (new): "+ str("%.2e" % U) )
  return Unc 

def energy_correction_near(rc, L, alpha):
  R = 40
  cor = 0.0
  for k0 in range(-R,R):
    for k1 in range(-R,R):
      for k2 in range(-R,R):
        Lk = L * sqrt(k0*k0 + k1*k1 + k2*k2)
        if Lk > rc:
          cor = cor + erfc(alpha * Lk) / (Lk)
  return 0.5 * cor

def energy_correction_far(kc, L, alpha):
  R = 2 * int(ceil(kc))
  cor = 0.0
  for k0 in range(-R,R):
    for k1 in range(-R,R):
      for k2 in range(-R,R):
        k = k0*k0 + k1*k1 + k2*k2
        if k > kc**2:
#         if min(k0, k1, k2) < -kc or max(k0, k1, k2) >= kc:
          cor = cor + exp(-(pi/alpha/L)**2.0 * k) / k
  return 1.0 / (2.0*pi*L) * cor

def kc2M(kc):
  return int( 2 * ceil(kc) + 2 )

def err_tot(eps_near, eps_far):
  if split_tol == 0.5:
    return eps_near + eps_far
  else:
    return sqrt(eps_near*eps_near + eps_far*eps_far)

def print_errors(kc):
  print("\t\tnear\t\tfar\t\tboth\t\tnumerical")
  eps_near = kolper_near_potential(rc, L, Q2, alpha)
  eps_far  = kolper_far_potential(kc, L, Q2, alpha)
  print("abs Potential:\t"+ str('%.2e' % eps_near) +"\t"+ str('%.2e' % eps_far) +"\t"+ str('%.2e' % err_tot(eps_near, eps_far)) +"\t"+ str('%.2e' % abs_error_potential))

  eps_near = kolper_near_energy(rc, L, Q2, alpha, N)
  eps_far  = kolper_far_energy(kc, L, Q2, alpha, N)
  print("abs Energy:\t"+ str('%.2e' % eps_near) +"\t"+ str('%.2e' % eps_far) +"\t"+ str('%.2e' % err_tot(eps_near, eps_far)) +"\t"+ str('%.2e' % abs_error_energy))

  eps_near = kolper_near_field(rc, L, Q2, alpha)
  eps_far  = kolper_far_field(kc, L, Q2, alpha)
  print("abs Field:\t"+ str('%.2e' % eps_near) +"\t"+ str('%.2e' % eps_far) +"\t"+ str('%.2e' % err_tot(eps_near, eps_far)) +"\t"+ str('%.2e' % abs_error_field))

  eps_near = kolper_near_force(rc, L, Q2, alpha, N)
  eps_far  = kolper_far_force(kc, L, Q2, alpha, N)
  print("abs Force:\t"+ str('%.2e' % eps_near) +"\t"+ str('%.2e' % eps_far) +"\t"+ str('%.2e' % err_tot(eps_near, eps_far)) +"\t"+ str('%.2e' % abs_error_force))

  eps_near = kolper_near_total_energy(rc, L, Q2, alpha)
  eps_far  = kolper_far_total_energy(kc, L, Q2, alpha)
  print("abs tot. Energy:"+ str('%.2e' % eps_near) +"\t"+ str('%.2e' % eps_far) +"\t"+ str('%.2e' % err_tot(eps_near, eps_far)) +"\t"+ str('%.2e' % abs_error_tenergy))

  print("\n\t\tnear\t\tfar\t\tboth\t\tnumerical")
  print("rel Potential:\t"+ str('%.2e' % 0) +"\t"+ str('%.2e' % 0) +"\t"+ str('%.2e' % 0) +"\t"+ str('%.2e' % rel_error_potential))
  print("rel Energy:\t"+ str('%.2e' % 0) +"\t"+ str('%.2e' % 0) +"\t"+ str('%.2e' % 0) +"\t"+ str('%.2e' % rel_error_energy))
  print("rel Field:\t"+ str('%.2e' % 0) +"\t"+ str('%.2e' % 0) +"\t"+ str('%.2e' % 0) +"\t"+ str('%.2e' % rel_error_field))
  print("rel Force:\t"+ str('%.2e' % 0) +"\t"+ str('%.2e' % 0) +"\t"+ str('%.2e' % 0) +"\t"+ str('%.2e' % rel_error_force))
  print("rel tot. Energy:"+ str('%.2e' % 0) +"\t"+ str('%.2e' % 0) +"\t"+ str('%.2e' % 0) +"\t"+ str('%.2e' % rel_error_tenergy))

def print_parameter_list():
  print("Options:")
  print("\tnd\t\t0d, 1d, 2d, or 3d")
  print("\ttestcase\t1 (cloud wall), 2 (silica melt), 3 (uniformly random), 4 (nacl crystal)")
  print("\ttestsize\tmultiplier of each box length")
  print("\ttol\t\taccuracy goal")
  print("\ttol-type\ttune for accuracy of 1 (abs. potential), 2 (abs. energy), 3 (abs. field), 4 (abs. force)")
  print("\tm\t\tNFFT window cutoff")
  print("\tndft\t\tchoose method for computing nonequispaced discrete Fourier transform:")
  print("\t\t\tndft (direct computation), nfft (fast approximation)")
  print("\tP\t\tnumber of gridpoints for regularization (ignored for 3dp)")
  print("\tp\t\tdegree of smoothness for regularization (ignored for 3dp)")
  print("\twindow\t\tNFFT window 1 (gaussian), 2 (bspline), 3 (sinc), 4 (kaiser), 5 (bessel_i0)")
  print("\tnon-periodic oversampling")
  print("\t\t\toversampling for non-periodic dimensions with M gridpoints:")
  print("\t\t\t-2 (double: 2*M), -1 (constant: M +2*m+2), value >= 0 (add: M + value)")
  print("\tperiodic oversampling")
  print("\t\t\toversampling for periodic dimensions with M gridpoints:")
  print("\t\t\t-2 (double: 2*M), -1 (constant: M +2*m+2), value >= 0 (add: M + value)")
  print("\tuse_kc\t\tuse radial Fourier space cutoff")
  print("\tignore_field\tdo not compute field")
  print("\tnprocs\t\tstart test runs with: mpirun -np nprocs\t")


if len(sys.argv) < 8:
  print("Usage: "+ os.path.basename(__file__) +" nd testcase testsize tol tol-type m ndft [P] [p] [window] [non-periodic oversampling] [periodic oversampling] [use_kc] [ignore_field] [nprocs]")
  print_parameter_list()
  sys.exit()

# We do not have a tuning for these parameters 
nd = set_str(1, "3dp")
tc = set_int(2, 3)
ts = set_int(3, 1)
tol = set_float(4, 5e-6)
tol_type = set_int(5, 1)

m  = set_int(6, 9)
ndft = set_str(7,"nfft")

if nd != "3d":
  if len(sys.argv) < 10:
    print("Usage: "+ os.path.basename(__file__) +" nd testcase testsize tol tol-type m ndft P p [window] [non-periodic oversampling] [periodic oversampling] [use_kc] [ignore_field] [nprocs]")
    print_parameter_list()
    sys.exit()

P = set_int(8, 18)
p = set_int(9, 9)

w = set_int(10, 2)
osn = set_number(11, 1)
osp = set_number(12, 0)
use_kc = set_str(13,"0")
ignore_field = set_str(14,"0")
nprocs = set_int(15,1)

# parameter sanity checks
if nd not in ("0d", "1d", "2d", "3d"):
  print("1st parameter 'nd' must be either '0d', '1d', '2d', or '3d'")
  sys.exit()
if tc not in (1,2,3,4):
  print("2nd parameter 'testcase' must be either '1', '2', '3', or '4'")
  sys.exit()
if ts < 1:
  print("3rd parameter 'testsize' must be at least 1")
  sys.exit()
if not tol > 0.0:
  print("4th parameter 'tol' must be positive")
  sys.exit()
if tol_type not in range(1,5):
  print("5th parameter 'tol-type' must be '1' (abs_pot), '2' (abs_energy), '3' (abs_field), or '4' (abs_force)")
  sys.exit()
if m < 1:
  print("6th parameter 'm' must be at least 1")
  sys.exit()
if ndft not in ("ndft","nfft"):
  print("7th parameter 'nfft' must be either \"ndft\", or \"nfft\"")
  sys.exit()
if P < 0:
  print("8th parameter 'P' must be non-negative")
  sys.exit()
if p < 0:
  print("9th parameter 'p' must be non-negative")
  sys.exit()
if w not in range(1,6):
  print("10th parameter 'window' must be either '1' (gaussian), '2' (bspline), '3' (sinc), '4' (kaiser), or '5' (bessel_i0)")
  sys.exit()
if osn < -2:
  print("11th parameter 'non-periodic oversampling' must be either 'val' > 0 (M+val), '-1' (M + 2*m+2), or '-2' (2*M)")
  sys.exit()
if osp < -2:
  print("12th parameter 'periodic oversampling' must be either 'val' > 0 (M+val), '-1' (M + 2*m+2), or '-2' (2*M)")
  sys.exit()
if nprocs < 1:
  print("15th parameter 'nprocs' must be > 0")
  sys.exit()

cwd = os.getcwd()

if w == 1:
  window = "gaussian"
elif w == 2:
  window = "bspline"
elif w == 3:
  window = "sinc"
elif w == 4:
  window = "kaiser"
elif w == 5:
  window = "bessel_i0"

test_dir = "nonperiodic" if nd == "0d" else nd +"-periodic"

# switch test case
if tc==1:
  # cloud_wall
  scale = 2**(ts-1)
  N  = 300 * scale**3
  L  = 10.0 * scale
  rc = 4
  Q2 = N
  Q4 = N
  tn = "cloud_wall"
  testname = "systems/"+ test_dir +"/cloud_wall_" + str(N) + ".xml.gz"
elif tc==2:
  # silica melt
  scale = 2**(ts-1)
  N  = 12960 * scale**3
  L  = 62.05966799605923 * scale
  rc = 8.0
  Q2 = 3.73248000e+04 * scale**3
  Q4 = 1.61243136e+05 * scale**3
  tn = "silica_melt"
  testname = "systems/"+ test_dir +"/silica_melt_" + str(N) + ".xml.gz"
elif tc==3:
  # random distribution
  N  = 11000 *(ts-1) + 1000
  L  = pow(N/1000.0, 1.0/3.0)
  rc = 0.62
#   rc = 1.2
  Q2 = N
  Q4 = N
  tn = "random_dist"
  testname = "numerics/accuracy-tests-increased-box/rand_"+ nd +"p_" + str(N) + "_fmm.xml.gz"
else:
  # nacl crystal
  N = 10648
  L = 22.0
  rc = 10.0
  Q2 = N
  Q4 = N
  tn = "nacl"
  testname = "systems/"+ test_dir +"/nacl_" + str(N) + ".xml.gz"

print("Tune for tol = "+ str('%.2e' % tol))
print("L  = "+ str('%.2f' % L))
print("rc = "+ str('%.2f' % rc))
print("Q2 = %.8e" % Q2)
# print("Q4 = %.8e" % Q4)

# invoke Kolafa/Perram parameter tuning
alpha = inv_kolper_alpha(rc, L, Q2, N, tol * split_tol, tol_type)
print("Kolafa-Perram Force = " + str(kolper_near_force(rc, L, Q2, alpha, N)))
kc    = inv_kolper_kc(alpha, L, Q2, N, tol * split_tol, tol_type)
M     = kc2M(kc)
print("alpha = "+ str('%.2f' % alpha))
print("kc = "+ str('%.2f' % kc))

# M = 64
# alpha = 0.7
# alpha = 0.748120

kc0 = kc
if nd == "2d":
  kc0 = inv_kolper_kc(alpha, 2*L, Q2, N, tol * split_tol, tol_type)
if nd == "1d":
  kc0 = inv_kolper_kc(alpha, sqrt(2)*2*L, Q2, N, tol * split_tol, tol_type)
if nd == "0d":
  kc0 = inv_kolper_kc(alpha, sqrt(3)*2*L, Q2, N, tol * split_tol, tol_type)

M0 = M1 = M2 = M
m0 = m1 = m2 = oversampled_gridsize(M, m, osp)
if nd == "2d":
  M0 = even( kc2M(kc0) + P )
  m0 = oversampled_gridsize(M0, m, osn)
if nd == "1d":
  M0 = M1 = even( kc2M(kc0) + P )
  m0 = m1 = oversampled_gridsize(M0, m, osn)
if nd == "0d":
  M0 = M1 = M2 = even( kc2M(kc0) + P )
  m0 = m1 = m2 = oversampled_gridsize(M0, m, osn)

h = L*M0/M

epsB = 0.5*P/M0
if nd == "3d":
  epsB = 0.0

# Compare correctness of near- and farfield error formulae.
# Therefore, choose either large rc or large kc.
# M = 2*M; kc = M/2.0-1.0
# rc = 1.5*rc

conf = " -c "
conf += "p2nfft_verbose_tuning,1,"
conf += "pfft_patience_name,estimate,"
conf += "p2nfft_ignore_tolerance,1,"
conf += "tolerance_field,1.000000e-20,"
conf += "p2nfft_r_cut,"+ str(rc) +","
conf += "p2nfft_alpha,"+ str(alpha) +","
conf += "pnfft_window_name,"+ window +","
conf += "pnfft_N,"+ str(M0) +","+ str(M1) +","+ str(M2) +","
conf += "pnfft_n,"+ str(m0) +","+ str(m1) +","+ str(m2) +","
conf += "pnfft_m,"+ str(m) +","
conf += "p2nfft_epsB,"+ str('%.4e' % epsB) +","
conf += "p2nfft_p,"+ str(p) +","
if ndft == "ndft":
  conf += "pnfft_direct,1,"
else:
  conf += "pnfft_direct,0,"

intpol_order = "3"
conf += "p2nfft_intpol_order,"+ intpol_order +","
conf += "pnfft_intpol_order,"+ intpol_order +","
conf += "pnfft_grad_ik,1,"
if use_kc != "0":
  conf += "p2nfft_k_cut,"+ str(kc) +","
if ignore_field != "0":
  conf += "p2nfft_ignore_field,1,"

conf += "p2nfft_reg_kernel,0,"

nofield=""
# nofield+=" -u nofield"

if nprocs > 1:
  startmpi = "mpirun -np "+ str(nprocs) +" "
else:
  startmpi = ""

line = startmpi + cwd + "/scafacos_test p2nfft " + testname + nofield + conf
print(line)

output = subprocess.check_output(line, stderr=subprocess.STDOUT, shell=True)
output = output.split('\n')
# output= ""

abs_error_potential = search_val(output, "abs_rms_potential_error", "=", "")
abs_error_energy    = search_val(output, "abs_rms_energy_error",    "=", "")
abs_error_field     = search_val(output, "abs_rms_field_error",     "=", "")
abs_error_force     = search_val(output, "abs_rms_force_error",     "=", "")
abs_error_tenergy   = search_val(output, "abs_total_energy_error",  "=", "")

rel_error_potential = search_val(output, "rel_rms_potential_error", "=", "")
rel_error_energy    = search_val(output, "rel_rms_energy_error",    "=", "")
rel_error_field     = search_val(output, "rel_rms_field_error",     "=", "")
rel_error_force     = search_val(output, "rel_rms_force_error",     "=", "")
rel_error_tenergy   = search_val(output, "rel_total_energy_error",  "=", "")

Q2 = search_val(output, "Q^2 =", "=", "")
Q4 = search_val(output, "Q^4 =", "=", "")
print("Test computed Q2 = %.8e" % Q2)
print("Test computed Q4 = %.8e" % Q4)

print("Error estimates with kc = "+ str('%.2f' % kc) +" predict following errors:")
print_errors(kc)

# print("Error estimates with kc = M/2 = "+ str('%.2f' % (M/2.0)) +" predict following errors:")
# print_errors(M/2.0)

print("Error estimates with kc = M/2-1 = "+ str('%.2f' % (M/2.0 - 1.0)) +" predict following errors:")
print_errors(M/2.0-1.0)

# print("Error estimates with kc = M/2-1 = "+ str('%.2e' % (M/2.0-1)) +" predict following errors:")
# print_errors(M/2.0-1.0)

rt_tune = search_val(output, "Time:", ":", "")
rt_all  = search_val(output, "Average time:", ":", "")
rt_near = search_val(output, "Near field computation takes", "takes", "s")
rt_far  = search_val(output, "Far field computation takes", "takes", "s")

print("\nRuntimes:")
print("all\t\tnear\t\tfar\t\ttune")
print(str('%.4e' % rt_all) +"\t"+  str('%.4e' % rt_near) +"\t"+ str('%.4e' % rt_far) +"\t"+ str('%.4e' % rt_tune))

outfile="parameters_"+ tn +"_"+ nd +"p_"+ str('%.2e'%tol) +"_"+ ndft
if use_kc != "0":
  outfile += "_kc"
outfile += ".txt"
if not os.path.isfile(outfile):
  with open(outfile, 'a') as f:
    f.write("# generated with:\n# \n" \
        "N          "+\
        "Q^2       "+\
        "L         "+\
        "rc        "+\
        "alpha     "+\
        "  M0  "+\
        "  M1  "+\
        "  M2  "+\
        "  m0  "+\
        "  m1  "+\
        "  m2  "+\
        "epsB      "+\
        "h         "+\
        " m  "+\
        "   P  "+\
        " p  "+\
        "abs_energy  "+\
        "abs_force   "+\
        "abs_tenergy  "+\
        "rel_energy  "+\
        "rel_force   "+\
        "rel_tenergy  "+\
        "rt-run    "+\
        "rt-near   "+\
        "rt-far    "+\
        "rt-tune   "+\
        "ndft  "+\
        "kc        "+\
        "nprocs  "+\
        "\n")

if use_kc != "0":
  kcstr = str('%.2e  ' % kc)
else:
  kcstr = "-1.0      "

with open(outfile, 'a') as f:
  f.write( \
      str('%9d  ' % N) +\
      str('%.2e  ' % Q2) +\
      str('%.2e  ' % L)  +\
      str('%f  ' % rc) +\
      str('%f  ' % alpha) +\
      str('%4d  ' % M0) +\
      str('%4d  ' % M1) +\
      str('%4d  ' % M2) +\
      str('%4d  ' % m0) +\
      str('%4d  ' % m1) +\
      str('%4d  ' % m2) +\
      str('%.2e  ' % epsB) +\
      str('%.2e  ' % h)    +\
      str('%2d  ' % m) +\
      str('%4d  ' % P) +\
      str('%2d  ' % p) +\
      str('%.2e    ' % abs_error_energy) +\
      str('%.2e    ' % abs_error_force)  +\
      str('%.2e     ' % abs_error_tenergy)+\
      str('%.2e    ' % rel_error_energy) +\
      str('%.2e    ' % rel_error_force)  +\
      str('%.2e     ' % rel_error_tenergy)+\
      str('%.2e  ' % rt_all)  +\
      str('%.2e  ' % rt_near) +\
      str('%.2e  ' % rt_far)  +\
      str('%.2e  ' % rt_tune) +\
      ndft + "  " +\
      kcstr +\
      str('%6d  ' % nprocs) +\
      "\n".expandtabs(2))

