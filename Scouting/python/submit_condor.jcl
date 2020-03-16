universe = vanilla
initialdir = .
use_x509userproxy = true
+C7OK = "yes"
Requirements = (C7 == "yes")
error =./con_logs/run_$(Process).error
log =./con_logs/run_$(Process).log
output =./con_logs/run_$(Process).out
executable = submit_condor.sh
arguments = $(Process)
Notification=never
queue 1000