#
# Sample QUERY_BC parsing script
#
# to run:
# > sh query_bc.sh 
#
# This shell script uses the linux command line utilities "grep", "paste", "awk", and "sed"
# to parse QUERY_BC lines from log files (grep), combine them into a single 
# table (paste), and then add the z-torque column in each row to produce the total
# torque vs time from all QUERY_BC zones (awk).
#
# note that we also use awk to zero the initial time (the first column) in each
# file...
#

#output log: rpm100wm/log.7172737.excalibur-wlm1.006...
grep "\[QUERY_BC:rotor:rotor]" rpm100wm/log.7172737.excalibur-wlm1.006 > rotor
grep "\[QUERY_BC:rotor:z1]" rpm100wm/log.7172737.excalibur-wlm1.006 > z1
grep "\[QUERY_BC:rotor:z0]" rpm100wm/log.7172737.excalibur-wlm1.006 > z0
grep "\[QUERY_BC:rotor:blade]" rpm100wm/log.7172737.excalibur-wlm1.006 > blade
grep "\[QUERY_BC:rotor:wall]" rpm100wm/log.7172737.excalibur-wlm1.006 > wall
paste rotor z1 z0 blade wall | awk '{print $3,$3-0.01372004,$9,$20,$31,$42,$53,$9+$20+$31+$42+$53}' > rpm100_torque.dat

#output log: rpm075wm/log.173.excalibur-wlm1.002...
grep "\[QUERY_BC:rotor:rotor]" rpm075wm/log.173.excalibur-wlm1.002 > rotor
grep "\[QUERY_BC:rotor:z1]" rpm075wm/log.173.excalibur-wlm1.002 > z1
grep "\[QUERY_BC:rotor:z0]" rpm075wm/log.173.excalibur-wlm1.002 > z0
grep "\[QUERY_BC:rotor:blade]" rpm075wm/log.173.excalibur-wlm1.002 > blade
grep "\[QUERY_BC:rotor:wall]" rpm075wm/log.173.excalibur-wlm1.002 > wall
paste rotor z1 z0 blade wall | awk '{print $3,$3-0.00212004,$9,$20,$31,$42,$53,$9+$20+$31+$42+$53}' > rpm075_torque.dat

#output log: rpm050wm/log.173.excalibur-wlm1.003...
grep "\[QUERY_BC:rotor:rotor]" rpm050wm/log.173.excalibur-wlm1.003 > rotor
grep "\[QUERY_BC:rotor:z1]" rpm050wm/log.173.excalibur-wlm1.003 > z1
grep "\[QUERY_BC:rotor:z0]" rpm050wm/log.173.excalibur-wlm1.003 > z0
grep "\[QUERY_BC:rotor:blade]" rpm050wm/log.173.excalibur-wlm1.003 > blade
grep "\[QUERY_BC:rotor:wall]" rpm050wm/log.173.excalibur-wlm1.003 > wall
paste rotor z1 z0 blade wall | awk '{print $3,$3-0.00228004,$9,$20,$31,$42,$53,$9+$20+$31+$42+$53}' > rpm050_torque.dat

#
# The first line of each of the 3 resulting files confirms the 
# columns being reported. You will need to go in afterwards and either remove or
# add a leading "#" to prevent this line from messing up certain plotting 
# programs.
#
# here is an easy way to add '# ' to the first line of each file using "sed":

sed -i '1s/^/# /' rpm100_torque.dat
sed -i '1s/^/# /' rpm075_torque.dat
sed -i '1s/^/# /' rpm050_torque.dat

