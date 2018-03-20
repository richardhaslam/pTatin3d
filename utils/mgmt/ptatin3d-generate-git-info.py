import os
import re

# Repo id
# git config --get remote.origin.url

# Revision id
# git log -1 --pretty="%H"

# Revision infomation
# git log -1 --pretty="(%h) : %an (%ae) | %ad" --date=iso

def GitFound_WriteInfoHeader():

    spacer = " "

    hf = open("ptatin_git_version_info.h",'w')
    hf.write("\n")
    hf.write("#ifndef __ptatin_git_version_info_h__\n")
    hf.write("#define __ptatin_git_version_info_h__\n")

    hf.write("\n")
    hf.write("\n")

    gitstr = "#undef PTATIN_VERSION_CNTR_REPO"
    hf.write(gitstr + "\n")
    gitstr = "#undef PTATIN_VERSION_CNTR_REVISION"
    hf.write(gitstr + "\n")
    gitstr = "#undef PTATIN_VERSION_CNTR_LOG"
    hf.write(gitstr + "\n")

    hf.write("\n")

    os.system("git status > tmpgit.txt")
    #magic_git_message "Your branch is up-to-date with"
    git_status = "\"git status: Your branch is up-to-date\""

    gitf = open("tmpgit.txt",'r')
    for line in gitf:
      if "Changes not staged for commit" in line:
        git_status = "\"git status: --- WARNING --- Your branch contains uncommitted changes\""
    gitf.close();

    gitstr = spacer.join( ["#define PTATIN_GIT_REPO_STATUS" , git_status ])
    hf.write(gitstr + "\n")

    hf.write("\n")
    # ---------------------------------------------------------------------
    os.system("git config --get remote.origin.url > tmpgit.txt")
    gitf = open("tmpgit.txt",'r')
    for line in gitf:
        giturl = line.rstrip('\n')
    gitf.close();

    giturl = re.sub("https://.*@","https://",giturl) # strip username from HTTPS URLs
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


def main():

    gitrepo_detected = False

    root = os.curdir
    for entry in os.listdir(root):
        if entry == ".git":
            gitrepo_detected = True

    if gitrepo_detected == True:
        GitFound_WriteInfoHeader()
    else:
        print('Unable to generate automated version header as git was not found')

if __name__ == "__main__":
    main()
