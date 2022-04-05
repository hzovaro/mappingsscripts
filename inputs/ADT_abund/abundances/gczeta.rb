#!/usr/bin/ruby
#
# v1.0.10
#
# Tool to make GC abundance files from a master file and a list of zetas
# Master file includes fiducial abundances plus Xi and  breaks on
# an Fe scale.   Currently in GC_Masterv2.6.3.txt
#
# Input files can be Linear or Log zetas and in any base zeta scaling,
# selected in the first two lines followed by a simle list of zetas.
# See sample zeta_xxx.txt files.
#
# ./gczeta.rb list
#
# eg
# ./gczeta.rb zeta_linearFe.txt
#
########################################################################
#
# setup and init data
#
gc_lines = File.readlines("GC_Masterv2.6.3.txt")
#
chi0 = -1.0
chi1 = +0.5
z_scale0 = "Fe"
#
z_list = ARGV[0]
zeta_lines = File.readlines("#{z_list}")
zeta_type = "Linear"   # or Log in input file header
#
n_elems = 30
el_num = {
  "H " => 1,  "He" => 2,  "Li" => 3,  "Be" => 4,   "B " => 5,
  "C " =>  6,  "N " => 7,  "O " => 8, "F " =>  9, "Ne" => 10,
  "Na" => 11, "Mg" => 12, "Al" => 13, "Si" => 14, "P " => 15,
  "S " => 16, "Cl" => 17, "Ar" => 18, "K " => 19, "Ca" => 20,
  "Sc" => 21, "Ti" => 22, "V " => 23, "Cr" => 24, "Mn" => 25,
  "Fe" => 26, "Co" => 27, "Ni" => 28, "Cu" => 29, "Zn" => 30 }
#
el_name = {
  1 => "H ",  2 => "He",  3 => "Li",  4 => "B ",  5 => "Be",
  6 => "C ",  7 => "N ",  8 => "O ",  9 => "F ", 10 => "Ne",
  11 => "Na", 12 => "Mg", 13 => "Al", 14 => "Si", 15 => "P " ,
  16 => "S ", 17 => "Cl", 18 => "Ar", 19 => "K ", 20 => "Ca",
  21 => "Sc", 22 => "Ti", 23 => "V ", 24 => "Cr", 25 => "Mn",
  26 => "Fe", 27 => "Co", 28 => "Ni", 29 => "Cu", 30 => "Zn"  }
#
# get breaks, Xi() and base abundances from master file
#
puts
#
xi  = Array.new(n_elems)
xi0 = Array.new(n_elems)
gc0 = Array.new(n_elems)
line_idx = 0
gc_lines.each do |l|
  line_idx += 1
  if line_idx <= 5
    puts " " << l
  end
  if line_idx == 3
    line = l.split(' ')
    z_scale0 = line[2]
    #    puts "#{z_scale0}"
  end
  if line_idx == 4
    line = l.split(' ')
    chi0 = Float(line[4])
    #    puts "#{chi0}"
  end
  if line_idx == 5
    line = l.split(' ')
    chi1 = Float(line[4])
    #    puts "#{chi1}"
  end
  if line_idx > 6
    line = l.split(' ')
    gc0[line_idx-7] = Float(line[2])
    el = Integer(line[0])
    xi0[line_idx-7] = 0.0
    if  el == 6
      xi0[line_idx-7] = Float(line[3])
    end
    if  el >= 8
      xi0[line_idx-7] = Float(line[3])
    end
  end
end
he_slope = 10.0**(gc0[el_num["He"]-1]+1.0783)-1.0
puts " He Slope: #{he_slope}"
#
# z_scale0 = element symbol for master file scale
# z_in_el  = element symbol for input file scale
#
#
# scale0_id   = base element idx of master file, usually Fe
# scale_in_id = base element idx of input file, usually Fe or O but may be any
# scaleO_id   = oxygen element idx for always haging zetaO for He
#
# zeta_in  = zeta value read in, log10 value always
# zeta_0   = zeta value in master file scale, usually Fe
# zeta_O   = zeta value in O scale always
#
scale0_id  = el_num[z_scale0] - 1  # scale used in Master data file
#
########################################################################
#  Log or Linear zeta inputs
#
zeta_type = zeta_lines[0].split(' ')[1]
########################################################################
# chosen element scale in input file
#
z_in_el = zeta_lines[1].split(' ')[1]
if (z_in_el.length == 1) then
  z_in_el += " "
end
#
scale_in_id = el_num[z_in_el] - 1  # 0 based for arrays
scaleO_id   = el_num["O "] - 1     # Oxygen scale independent of input
scaleFe_id  = el_num["Fe"] - 1    # Iron scale independent of input
#
# get master Xi for chosen input scale element
#
xi_in   = xi0[scale_in_id]
puts " Xi(#{z_in_el})_#{z_scale0} = #{xi_in}"
puts
puts "   ========================"
#
# get scaled cuts from master => input scale
#
chi_lo  = chi0 + xi0[scale_in_id]
chi_hi  = chi_lo * ( chi1/chi0 )
#
puts "      Chi_lo_#{z_in_el} = " << format("%6.3f", chi_lo)
puts "      Chi_hi_#{z_in_el} = " << format("%6.3f", chi_hi)
puts "      Xi(#{z_scale0})_#{z_in_el} = " << format("%6.3f", xi0[scale0_id] - xi0[scale_in_id])
#
puts "   ========================"
#
# change xi to new base
#
specials = [1, 2, 3, 4, 5, 7] # H, He, Li, Be, B, N
n_elems.times do |idx|
  xi[idx] = xi0[idx] - xi_in
  if specials.include?(idx+1) then
    xi[idx] = 0.0
  end
#  puts "      Xi(#{el_name[idx+1]})_#{z_in_el} = "<< format("%6.3f",xi[idx])
end
#
#puts "   ========================"
puts
########################################################################
#
# delta function for scale conversion and abundances offset
#
def delta_zxicc( z, x, id, c_lo, c_hi )
  del = 0.0
  if ( z > c_hi ) then
    del = ( x[ id ]/c_lo )*c_hi
  elsif ( z > c_lo ) then
    del = ( x[ id ]/c_lo )*z
  else
    del = x[ id ]
  end
  return del
end
#
########################################################################
#
# get zetas to evaluate and write to stdout
#
line_idx  = 0
#
zeta_in   = 0.0  # from input file line
zeta_0    = 0.0  # in master file scale
zeta_O    = 0.0  # in oxygen scale
linzeta_O = 0.0
#
zeta_lines.each do |l|
  #
  line_idx += 1
  if ( line_idx > 2 ) then # skip two line this time..
    ########################################################################
    #
    # read in zeta_in, also make zeta_0 (main GC data file scale) and
    # zeta_O for helium, independent of zeta_in scale.
    #
    # all zetas in dex unless otherwise noted
    #
    puts
    #
    #
    if (zeta_type == "Linear") then
      zeta_in = Math.log10(Float(l))  # read in a linear zeta, cvt to log
    else
      zeta_in = Float(l)  # read in a log zeta
    end
    linzeta_in = 10.0**zeta_in
    ########################################################################
    # open an output abn file
    #
    file_prefix = "GC_Z" << (z_in_el).strip << "_"
    file_name = "GCZ_XXX.abn"
    jfile_name = "GCZ_MXXX.json"
    if (zeta_type == "Linear") then
      file_id = Integer((linzeta_in*1000).round)
      file_name  = file_prefix + format("%4.4d",file_id) + ".abn"
      jfile_name = file_prefix + format("%4.4d",file_id) + ".json"
    else
      file_id = Integer((zeta_in*1000).round)
      if (file_id < 0) then
        file_name  = file_prefix + "M" + format("%4.4d",-file_id) + ".abn"
        jfile_name = file_prefix + "M" + format("%4.4d",-file_id)  + ".json"
      else
        file_name = file_prefix + "P" + format("%4.4d",file_id) + ".abn"
        jfile_name = file_prefix + "P" + format("%4.4d",file_id) + ".json"
      end
    end
    puts file_name
    # ADT: I don't want the json file, so commenting this out#   puts jfile_name
    #
    # create abund and JSON file at the same time
    #
    File.open(file_name, 'w') do |f|
    #
     # ADT: File.open(jfile_name, 'w') do |fj|
     # ADT: fj.puts "{"

      f.puts "\%"
      f.puts "\% GC2016 Zeta-Delta Scaled abundance file"
      f.puts "\% Generated by gczeta v1.0.10"
      dt    = `date`
      f.puts "\% Generated  :#{dt}"
      f.puts "\%"
      f.puts "\% Master File: GC_Masterv2.6.3.txt"
      f.puts "\% Master Base Element: " << z_scale0
      f.puts "\%     Low Cut: " << format("%6.3f", chi0) <<" dex"
      f.puts "\%    High Cut: " << format("%6.3f", chi1) <<" dex"
      f.puts "\%"
      f.puts "\%   File Base Element: " << z_in_el
      f.puts "\%     Low Cut: " << format("%6.3f", chi_lo) <<" dex"
      f.puts "\%    High Cut: " << format("%6.3f", chi_hi) <<" dex"
    #
    ########################################################################
    #
    # make other zetas:
    #
    delta = delta_zxicc( zeta_in, xi, scale0_id, chi_lo, chi_hi )
    zeta_0 = delta + zeta_in # in master scale
    #
    delta = delta_zxicc( zeta_in, xi, scaleO_id, chi_lo, chi_hi )
    scaled_O  = gc0[scaleO_id] + delta + zeta_in
    zeta_O    = delta + zeta_in # oxygen scale
    #
    delta = delta_zxicc( zeta_in, xi, scaleFe_id, chi_lo, chi_hi )
    scaled_Fe  = gc0[scaleFe_id] + delta + zeta_in
    zeta_Fe   = delta + zeta_in # iron scale
    #
    linzeta_O = 10.0**zeta_O  # for helium scaling
    scaled_He = -1.0783+Math.log10(1.0+he_slope*linzeta_O)

    f.puts "\% "
    f.puts "\% File abundance: log zeta_" << el_name[scale_in_id+1] << " = " << format("%7.4f",zeta_in) <<" dex"
    f.puts "\%               : linear zeta_" << el_name[scale_in_id+1] << " = " << format("%7.4f",linzeta_in)
    f.puts "\% "
    f.puts "\% Absolute Abundances:"
    f.puts "\%          log(He/H) = " << format("%7.4f",scaled_He) <<" dex"
    f.puts "\%          log(O/H)  = " << format("%7.4f",scaled_O) <<" dex"
    f.puts "\%          log(Fe/H) = " << format("%7.4f",scaled_Fe) <<" dex"
    f.puts "\% "
    f.puts "\%       12+log(He/H) = " << format("%7.4f",12.0 + scaled_He) <<" dex"
    f.puts "\%       12+log(O/H)  = " << format("%7.4f",12.0 + scaled_O) <<" dex"
    f.puts "\%       12+log(Fe/H) = " << format("%7.4f",12.0 + scaled_Fe) <<" dex"
    f.puts "\% "
    f.puts "\% GC2016 Fe :   log zeta_Fe = " << format("%7.4f",zeta_Fe) <<" dex"
    f.puts "\% GC2016 O  :    log zeta_O = " << format("%7.4f",zeta_O ) <<" dex"
    f.puts "\% "
    f.puts "\% Comparison relative abundances at the same absolute abundances:"
    f.puts "\% "
    f.puts "\% AGSS09 Fe : AGSS09 [Fe/H] = " << format("%7.4f",zeta_Fe + (4.50 - 4.48)) <<" dex"
    f.puts "\% AG1989 Fe : AG1989 [Fe/H] = " << format("%7.4f",zeta_Fe + (4.33 - 4.48)) <<" dex"
    f.puts "\% "
    f.puts "\% AGSS09 O  : AGSS09  [O/H] = " << format("%7.4f",zeta_O  + (3.31 - 3.24)) <<" dex"
    f.puts "\% AG1989 O  : AG1989  [O/H] = " << format("%7.4f",zeta_O  + (3.07 - 3.24)) <<" dex"
    f.puts "\% "
    if (zeta_type == "Linear") then
      f.puts " GC2016 Zeta Scaled abundances, linear zeta_" << el_name[scale_in_id+1] << " = " << format("%7.4f",linzeta_in)
    else
      f.puts " GC2016 Zeta Scaled  abundances, log zeta_" << el_name[scale_in_id+1] << " = " << format("%7.4f",zeta_in)
    end
    f.puts " #{n_elems} 1"
    #
    puts "========================================"
    puts " log zeta_#{z_scale0}  log zeta_#{z_in_el}  log zeta_O "
    puts "========================================"
    puts format("   %9.6f",zeta_0) << format("    %9.6f",zeta_in) << format("   %9.6f",zeta_O)
    ########################################################################
    #
    # now scale all the elements with zeta_in
    #
    puts "========================================"
    puts " El  log[GC]   Xi(El)    Delta   logNew"
    puts "========================================"
    #
    # create JSON file for website at the same time
    #

    n_elems.times do |i|
      #
      delta  = delta_zxicc( zeta_in, xi, i, chi_lo, chi_hi )
      gc_new = gc0[i] + delta + zeta_in  # input file scale
      #
      # now special case some elements
      #
      if specials.include?(i+1) then
        if (i == 0) then # Hydrogen = 0.0
          gc_new = 0.0
        end
        if (i == (el_num["He"]-1)) then # helium scales linearly with zeta_O
          gc_new = scaled_He
        end
        if (i == (el_num["Li"]-1)) then # Li constant
          gc_new = gc0[2]
        end
        if (i == (el_num["Be"]-1)) then # Be constant
          gc_new = gc0[3]
        end
        if (i == (el_num["B "]-1)) then # B constant
          gc_new = gc0[4]
        end
        if (i == (el_num["N "]-1)) then # Nitrogen scaling
          #gc_new = scaled_O +Math.log10(10.0**(-1.732) + 10.0**(2.19+scaled_O))
          delta = Math.log10(10.0**(-0.764) + 10.0**(zeta_O-0.082))
          gc_new = gc0[i] + delta + zeta_O
        end
      end
      #
      # write out, tidying up longer numbers when <= -10.000
      #
      if ( gc0[i] < -10.0) then
        if ( gc_new < -10.0) then
          puts " #{el_name[i+1]} "<< format("  %5.2f",gc0[i]) << format("  %7.4f",xi[i]) << format("  %7.4f",delta) << format("  %6.3f",gc_new)
        else
          puts " #{el_name[i+1]} "<< format("  %5.2f",gc0[i]) << format("  %7.4f",xi[i]) << format("  %7.4f",delta) << format("  %7.4f",gc_new)
        end
        else
        if ( gc_new < -10.0) then
          puts " #{el_name[i+1]} "<< format("  %6.3f",gc0[i]) << format("  %7.4f",xi[i]) << format("  %7.4f",delta) << format("  %6.3f",gc_new)
        else
          puts " #{el_name[i+1]} "<< format("  %6.3f",gc0[i]) << format("  %7.4f",xi[i]) << format("  %7.4f",delta) << format("  %7.4f",gc_new)
        end
      end
      #
      if ( gc_new < -10.0) then
        f.puts format(" %2d",(i+1)) << format("  %6.3f",gc_new)
        else
        f.puts format(" %2d",(i+1)) << format("  %7.4f",gc_new)
      end
      #
       # ADT: if ( gc_new < -10.0) then
       # ADT:   fj.puts "\"#{el_name[i+1]}\"" << format(":  %6.3f,",gc_new)
       # ADT:   else
       # ADT:   fj.puts "\"#{el_name[i+1]}\"" << format(":  %7.4f,",gc_new)
       # ADT: end
      #
    end
    puts "========================================"
     # ADT: fj.puts "}"
     # ADT: end # output jfile
    end # output file
  end
end
puts


