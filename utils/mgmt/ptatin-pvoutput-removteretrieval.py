
import os
import sys
import getopt

def usage():
	print '[ptatin-pvoutput-remoteretrieval]: scp tool to extract ParaView output files and pTat3d logfile'
	print 'Usage:'
	print '  -t TIMESTEP (or --timestep=TIMESTEP) : specify the output file you want'
	print '  -m MACHINE  (or --machine=MACHINE)   : specify the machine where the files live'
	print '  -p PATH     (or --path=PATH)         : specify the path to the files'


def main():
	print '-- ptatin-pvoutput-remoteretrieval --'

	# parse options
	try:
		opts, args = getopt.getopt(sys.argv[1:], "ht:p:m:", ["timestep=", "path=", "machine="])
	except getopt.GetoptError as err:
		print 'Error: Unknown option provided'
		usage()
		sys.exit(2)

	if len(opts) == 0:
		print 'Error: No options provided'
		usage()
		sys.exit(2)

	timestep = 0
	path = 0
	machine = "dmay@musashi.ethz.ch"

	for opt, arg in opts:
		if opt in ("-t", "--timestep"):
			timestep = arg
		elif opt in ("-p", "--path"):
			path = arg
		elif opt in ("-m", "--machine"):
			machine = arg
		elif opt in ("-h"):
			usage()


	print '  Fetching timestep:', timestep
	print '  From:', machine, '-->>', path


	files_to_fetch =                  path + "/" + "step" + timestep + "* "
	files_to_fetch = files_to_fetch + path + "/*.pvts "
	files_to_fetch = files_to_fetch + path + "/*.pvtu "
	files_to_fetch = files_to_fetch + path + "/*.pvd "
	files_to_fetch = files_to_fetch + path + "/ptatin.options "
	files_to_fetch = files_to_fetch + path + "/ptatin.log* "
	files_to_fetch = files_to_fetch + path + "/ptatin.petsc.log_summary* "
	#print files_to_fetch

	cmd = 'scp ' + machine + ':\"' + files_to_fetch + '\" .'
	#print cmd
	os.system(cmd)


if __name__ == "__main__":
	main()
