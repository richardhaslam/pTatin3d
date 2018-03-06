
import os

# Repo id
# git config --get remote.origin.url

# Revision id
# git log -1 --pretty="%H"

# Revision infomation
# git log -1 --pretty="(%h) : %an (%ae) | %ad" --date=iso

def pTatin3d_GitFound_WriteInfoHeader():

    spacer = " "
    
    hf = open("ptatin_version_info.h",'w')
    hf.write("\n")
    hf.write("#ifndef __ptatin_version_info_h__\n")
    hf.write("#define __ptatin_version_info_h__\n")

    hf.write("\n")
    hf.write("#define PTATIN_DEVELOPMENT_VERSION\n")
    hf.write("\n")
    hf.write("\n")

    # ---------------------------------------------------------------------
    os.system("git config --get remote.origin.url > tmpgit.txt")
    gitf = open("tmpgit.txt",'r')
    for line in gitf:
        giturl = line.rstrip('\n')
    gitf.close();

    # We do this splitting nonsense to remove the username from the URL
    # Note that this will only work with HTTPS addresses (not SSH ones)
    tmp = giturl.split("bitbucket.org")
    giturl = "https://bitbucket.org" + tmp[1]
    giturl = '"%s"' % ("git url: " + giturl)
    gitstr = spacer.join( ["#define PTATIN_VERSION_CNTR_REPO" , giturl ])
    hf.write(gitstr + "\n")

    # ---------------------------------------------------------------------
    os.system("git log -1 --pretty=\"%H\" > tmpgit.txt")
    gitf = open("tmpgit.txt",'r')
    for line in gitf:
        gitrev = line.rstrip('\n')
    gitf.close();

    gitbranch = []
    os.system("git branch > tmpgit.txt")
    gitf = open("tmpgit.txt",'r')
    for line in gitf:
      if line[0] == "*":
        gitbranch = line.replace('*',' ')
        gitbranch = gitbranch.lstrip()
        gitbranch = gitbranch.rstrip('\n')
    gitf.close();

    gitrev = '"%s"' % ("commit hash: " + gitrev + " (" + gitbranch + ")")
    gitstr = spacer.join( ["#define PTATIN_VERSION_CNTR_REVISION" , gitrev ])
    hf.write(gitstr + "\n")

    # ---------------------------------------------------------------------
    os.system("git log -1 --pretty=\"log[%h]: %an (%ae) | %ad\" --date=iso > tmpgit.txt")
    gitf = open("tmpgit.txt",'r')
    for line in gitf:
        gitlog = line.rstrip('\n')
    gitf.close();
    gitlog = '"%s"' % (gitlog)
    gitstr = spacer.join( ["#define PTATIN_VERSION_CNTR_LOG" , gitlog ])
    hf.write(gitstr + "\n")

    hf.write("\n")
    hf.write("#endif\n")
    hf.write("\n")

    hf.close()

    os.system("rm -f tmpgit.txt")


def pTatin3d_Default_WriteInfoHeader():
    
    spacer = " "
    
    hf = open("ptatin_version_info_default.h",'w')
    hf.write("\n")
    hf.write("#ifndef __ptatin_version_info_h__\n")
    hf.write("#define __ptatin_version_info_h__\n")
    
    hf.write("\n")
    hf.write("#define PTATIN_DEVELOPMENT_VERSION\n")
    hf.write("\n")
    hf.write("\n")
    
    # ---------------------------------------------------------------------
    gitstr = "#define PTATIN_VERSION_CNTR_REPO \"git url: https://bitbucket.org/jedbrown/ptatin3d.git\""
    hf.write(gitstr + "\n")
    
    # ---------------------------------------------------------------------
    gitstr = "#define PTATIN_VERSION_CNTR_REVISION \"commit hash: [out-of-date] Execute " + '\\"make releaseinfo\\"' + " to update to the most recent revision\""
    hf.write(gitstr + "\n")
    
    # ---------------------------------------------------------------------
    gitstr = "#define PTATIN_VERSION_CNTR_LOG \"log: [out-of-date] Execute " + '\\"make releaseinfo\\"' + " to update to the most recent revision\""
    hf.write(gitstr + "\n")
    
    hf.write("\n")
    hf.write("#endif\n")
    hf.write("\n")
    
    hf.close()

def pTatin3d_WriteVersionInfoHeader(hash,major,minor,patch):
    
    spacer = " "
    
    hf = open("ptatin_version_info.h",'w')
    hf.write("\n")
    hf.write("#ifndef __ptatin_version_info_h__\n")
    hf.write("#define __ptatin_version_info_h__\n")
    
    hf.write("\n")
    hf.write("#define PTATIN_RELEASE\n")
    hf.write("\n")
    hf.write("\n")
    
    # ---------------------------------------------------------------------
    gitstr = "#define PTATIN_VERSION_CNTR_REPO \"git url: https://bitbucket.org/jedbrown/ptatin3d.git\""
    hf.write(gitstr + "\n")
    
#    gitstr = "#define PTATIN_VERSION_CNTR_REVISION \"" + str(hash) + "\""
#    hf.write(gitstr + "\n")
    gitrev = '"%s"' % ("commit hash: " + hash )
    gitstr = spacer.join( ["#define PTATIN_VERSION_CNTR_REVISION" , gitrev ])
    hf.write(gitstr + "\n")

    # ---------------------------------------------------------------------
    gitstr = "#define PTATIN_VERSION_MAJOR " + str(major)
    hf.write(gitstr + "\n")
    gitstr = "#define PTATIN_VERSION_MINOR " + str(minor)
    hf.write(gitstr + "\n")
    gitstr = "#define PTATIN_VERSION_PATCH " + str(patch)
    hf.write(gitstr + "\n")
    
    hf.write("\n")
    hf.write("#endif\n")
    hf.write("\n")
    
    hf.close()


def main():
    # Default file which lives in the repo
    #pTatin3d_Default_WriteInfoHeader()
    #return	    

    gitrepo_detected = False
    
    root = os.curdir
    for entry in os.listdir(root):
        if entry == ".git":
            gitrepo_detected = True
    
    if gitrepo_detected == True:
        pTatin3d_GitFound_WriteInfoHeader()
    else:
        # Hard coded revision key and release tags
        git_commit_hash_key = "abcd01234"
        ptat_rev_major = 1
        ptat_rev_minor = 0
        ptat_rev_patch = 0
        
        pTatin3d_WriteVersionInfoHeader(git_commit_hash_key,ptat_rev_major,ptat_rev_minor,ptat_rev_patch)


if __name__ == "__main__":
    main()
