restart:
# Duffing Equation: diff( x, t$2 ) = -a x^3 + b (x-xasym) - c diff(x,t) + F0 cos(wt)

# Parameters
#################################
#                               #
#       Input parameters        #
#                               #
#################################

#   for one instance of CHAOS, set a=b=c=1, xasym=0, w=1, F0=0.900, x0=0, v0=0
#   for one instance of 3RD SUBHARMONIC, set a=b=c=1, xasym=0, w=1, F0=0.850, x0=0, v0=0

a := 1.000:        # coefficient of x^3 in duffing equation
b := 1.000:        # coefficient of x in duffing equation
c := 1.000:        # dimensionless damping parameter ( c >= 0 )
xasym := 0:        # potential asymmetric offset 

w  := 1.000:       # ratio of driving freq over natural frequency
F0 := 0.855:       # dimensional driving amplitude

x0 := 0.000:       #  initial (time = 0) angle (in terms of theta0 )
v0 := 0.000:       #  initial (time = 0)  angular speed  (1/sec)

################################
#     sampling parameters      #
################################

time_step := evalf(2*Pi/wd)/1: # time between samples, 2*Pi/wd for Poincare section
time_lag := 100*time_step:     # time at which to start sampling, time_step/4 for 90 degree lag
NPP := 1000:                   # number of poincare points
# Solve
# this takes a couple of minutes for 4000 points at Digits=12 precision

###############################################################
#                                                             #
#   Maple code which produces an animation of the             #
#   periodically-driven, damped, inverted pendulum.           #
#   The parameters of the problem are chosen in such a        #
#   way that there are two stable static equilibrium          #
#   angles  +theta0 and -theta0  with respect to the          #
#   vertical, and theta0 is sufficiently small such that      #
#           sin(theta0) = theta0 - theta0^3/6                 #
#   to good approximation. The variable "x" describes         #
#   the angle of the pendulum relative to theta0:             #
#           x = theta / theta0                                #
#   so +1 and -1 are the static equilibrium positions         #
#   for x.                                                    #
#                                                             #
###############################################################

###############################################################
#
#  The user need not change anything below this line.
#
###############################################################

Digits:=12:        # precision of calculations, changing this may be a bad idea

datascale := 800/3.:           # to scale data to resemble experimental output

Fd := F0:  wd := w:
c_  := .15*c:    # hack so that setting c := 1 above gives chaos
Fd_ := .30*Fd:     # hack so that setting Fd := 1 above gives chaos
  
time_end := time_lag+NPP*time_step:

#  perform some error checks

flag:=true:

if (time_step <= 0.0) then
   print("  ERROR: time_step > 0 required"):
   flag:=false:
fi:

if (c_ < 0) then
   print("  ERROR: dampening factor c >= 0 required"):
   flag:=false:
fi:

if (time_lag < 0.0) then
   print("  ERROR: time_lag >= 0 required"):
   flag:=false:
fi:

if (time_end <= time_lag ) then
   print("  ERROR: time_end > time_lag required"):
   flag:=false:
fi:

if (c_>0 and Fd_>0) then
   poincare_section_period:=evalf(2*Pi/wd):
else
   poincare_section_flag:=false:
fi:

if (flag) then

                 ##########################################
                 #                                        #
                 #     SOLVE THE EQUATIONS OF MOTION      #
                 #                                        #
                 ##########################################

##############################
#
#  the functions to determine  
#  (angle, angular speed)
#
##############################

funcs:={ x(t), v(t)  }:

############################
#
#  the equations of motion
#
##############################

sysODE:={
  D(x)(t) = v(t)  ,
  D(v)(t) = -a*x(t)^3 + b*(x(t)-xasym) - c_*v(t) + Fd_*cos(wd*t)
}:

##############################
#
#  initial conditions
#
##############################

ICs:={
  x(0) = x0,
  v(0) = v0
}:

##############################
#
#  Now solve the equations of motion!! First, set up the array
#  "tvals" which contains the values of time at which to compute
#  the functions.  "dsolve" returns the results in a specific
#  matrix, so then extract the results into a more useful form.
#  Put the values of x(t) into x_vals, etc.
#
##############################

Npoints:=round((time_end-time_lag)/time_step)+1;

print( `Digits `, Digits, ` Npoints`, Npoints );
tvals:=array(1..Npoints):
for i from 1 to Npoints do
   tvals[i]:=time_lag+(i-1)*time_step:
   od:

answer:=dsolve(sysODE union ICs, funcs, numeric, method=dverk78, output=tvals):

funcnames:=[seq(op(0,answer[1,1][i]),i=2..3)]:
member('x',funcnames,'ind'): ind:=ind+1:
x_vals:=array(1..Npoints):
for i from 1 to Npoints do
   x_vals[i]:=datascale*answer[2,1][i,ind]:
   od:
member('v',funcnames,'ind'): ind:=ind+1:
v_vals:=array(1..Npoints):
for i from 1 to Npoints do
   v_vals[i]:=datascale*answer[2,1][i,ind]:
   od:

#save x,v data to file.  multiply by datascale and round to get expt output.
fd := fopen("xduf.txt", WRITE):
for i from 1 to Npoints do: fprintf( fd, "%g\n ", round(x_vals[i]) ): od:
fclose(fd):
fd := fopen("vduf.txt", WRITE):
for i from 1 to Npoints do: fprintf( fd, "%g\n ", round(v_vals[i]) ): od: 
fclose(fd):

print(` Computation done `):
fi:
# Plot
# this also takes a couple of minutes for 4000 points at Digits=12 precision

         ##########################################################
         #                                                        # 
         # Plot the solution of the duffing differential equation #
         #                                                        #
         ##########################################################

# Global variables
#              x_vals,v_vals: Result points
# Local variables
#              x1 : list of position values of the points to plot
#              y1 : list of velocity values of the points to plot
#              pts : points in plotable format
x1 := convert(x_vals,list):
y1 := convert(v_vals,list):
pts := zip((x,y) -> [x,y], x1,y1):
llimit := 1.5*datascale:
plot(pts,-llimit..llimit,-llimit..llimit, style=point,symbol=circle,color=red);

# 
# Animate - UNDER CONSTRUCTION
#### !!!! UNDER CONSTRUCTION !!!! ####



                 ##########################################
                 #                                        #
                 #         PERFORM THE ANIMATION          #
                 #                                        #
                 ##########################################

print(` Preparing the animation `, Npoints):

   #######################################
   #
   #   Given angle "x", this draws
   #   a pendulum of unit length hanging
   #   at an angle of x wrt vertical.
   #   color_flag=0  --> blue pendulum
   #   color_flag<>0 --> yellow pendulum (useful for Poincare sections)
   #
   #######################################
 
draw_pendulum:=proc(x,color_flag)
 global colors_to_cycle,arrow_body_width,arrow_head_width,
        arrow_head_height,theta0:
 local origin,tip:
 origin:=[0,0]:
 tip:=[sin(x*theta0),cos(x*theta0)]:
 if (color_flag=0) then
    return PLOT(arrow(origin,tip,arrow_body_width,arrow_head_width,arrow_head_height,color=blue)):
 else
    return PLOT(arrow(origin,tip,arrow_body_width,arrow_head_width,arrow_head_height,color=colors_to_cycle[color_flag])):
 fi:
end:

   #  needed quantities to perform the animation

fscale:=1.0:
arrow_body_width:=0.05*fscale:
arrow_head_width:=0.0:
arrow_head_height:=0.0:
space_range:=1:
color_count:=0:
color_change_flag:=true:
colors_to_cycle:=array(1..8):
colors_to_cycle[1]:=yellow:
colors_to_cycle[2]:=green:
colors_to_cycle[3]:=red:
colors_to_cycle[4]:=cyan:
colors_to_cycle[5]:=magenta:
colors_to_cycle[6]:=navy:
colors_to_cycle[7]:=aquamarine:
colors_to_cycle[8]:=violet:

with(plottools):

cs:=proc(x)
 if (x=0) then "0"
 elif (abs(x)<1) then
    if (x>0) then
       cat("0",convert(x,string)):
    else
       cat("-0",convert(abs(x),string)):
    fi:
 else
    convert(x,string):
 fi:
end:

  #  routine used for the Poincare section strobing

cflag:=proc(tval)
 global poincare_section_period,time_lag,color_change_flag,
        time_step,color_count,number_of_points_in_section:
 if (abs(frem(tval-time_start_plots,poincare_section_period))<time_step) then
    if (color_change_flag) then
       color_count:=(color_count mod number_of_points_in_section)+1:
       color_change_flag:=false:
    fi:
    return 1;
 else
    color_change_flag:=true:
    return 0;
 fi:
end:


Pstring:=cat("c = ",cs(c),", wd = ",cs(wd),", R = ",cs(R)):
Istring:=cat("x(0) = ",cs(x0),",  v(0) = ",cs(v0)):
titlestring:=cat("Forced-damped inverted pendulum\n",Pstring,"\n",Istring):

plotxv := true:
if (not plotxv) then
if (poincare_section_flag) then
   frames:=[seq(draw_pendulum(x_vals[i],cflag(tvals[i])),i=1..Npoints)]:
else
   frames:=[seq(draw_pendulum(x_vals[i],0),i=1..Npoints)]:
fi:
xeqb:=1.3*sin(theta0):
yeqb:=1.3*cos(theta0):
animation:=plots[display](PLOT(CURVES([[-xeqb,yeqb],[0,0],[xeqb,yeqb]]),
                  COLOR(RGB,1,0,0)),
                  plots[display](frames,insequence=true),
                  axes=normal,tickmarks=[0,0],
                  title=titlestring,titlefont=[HELVETICA,16],
                  view=[-space_range..space_range,-space_range..space_range]):

else
if (poincare_section_flag) then
  frames:=[seq(PLOT(POINTS([x_vals[i],v_vals[i]],SYMBOL(CIRCLE)),COLOR(RGB,cflag(tvals[i]),0,1)),i=1..Npoints)]:
else
frames:=[seq(PLOT(POINTS([x_vals[i],v_vals[i]],SYMBOL(CIRCLE)),COLOR(RGB,0,0,1)),i=1..Npoints)]:
fi:
space_range := 5:
#animation:=plots[display]( PLOT(
                                CURVES( [[1,0],[0,0],[0,1]] ),
                                COLOR(RGB,1,0,0)
                               ),
                  plots[display](frames,insequence=true),
                  axes=normal,tickmarks=[0,0],
                  title=titlestring,titlefont=[HELVETICA,16],
                  view=[-space_range..space_range,-space_range..space_range]):
fi:
#print(animation):
# click on animation for play/stop/... controls at top of window

