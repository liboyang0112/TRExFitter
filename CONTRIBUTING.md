# Just show me the commands

For a quick svn-git translation table, see the
[git-svn crash course](https://git.wiki.kernel.org/index.php/GitSvnCrashCourse)

# How to contribute

This doc is mainly for developers, hopefully serving as a quick reference for using git.

GitLab repos can be used in numerous ways, but the recommended way is to set up your GitLab account with ssh keys and use the ssh protocol.

### Setup git

On lxplus, git is already installed. Just check what version you are using:
```bash
git --version
```

For first time use, you must setup your configuration,
[cern gitlab "getting started"](https://cern.service-now.com/service-portal/article.do?n=KB0003137)
In addition to username and email might be interested in setting your
preferred editor and colored diff. However, don't bother with the
password caching.  Instead, [setup your ssh
keys](https://gitlab.cern.ch/help/ssh/README.md) or use kerberos.

```bash
git config --global user.name "First Last"
git config --global user.email "first.last@cern.ch"
git config --global color.ui auto
git config --global core.editor "emacs -nw"

```

You might have to require us to add you as a developer before you can
push changes to the gitlab repository.  (but, in the meantime, you can
still commit them to your local repository)

### Setup a clean development area
```bash
mkdir TRexFitter
cd TRexFitter
git clone ssh://git@gitlab.cern.ch:7999/TRExStats/TRExFitter.git
```

### Commit changes to your local repository and push them to github:
```bash
cd TRexFitter
git status # start from a clean state, up-to-date with the origin
git pull   # make sure you are up to date
git checkout -b some-project  # put your work on a branch
git commit some.txt file.txt  # commit files and prompt editor for commit message
git push origin some-project  # push your changes to github
```
Notes:
- feature development should generally always occur on a dedicated
  branch, rather than on the master branch
- try to pick a branch name that is short and describes what you are
  trying to do
- when you have made enough progress, push to github, then issue a
  pull request (green button), and start discussing the proposed
  changes
- recommendations about commit messages. Follow Linus' guidelines: one
  short descriptive line, followed by a longer description of the
  motivation and implementation, if necessary. Limit line length (with
  emacs `M-q` or `M-x fill-paragraph`)
- setting up your git preferences (editor, aliases, etc.): see the
  [Git Configuration](http://git-scm.com/book/en/Customizing-Git-Git-Configuration)
  docs
- report issues and bugs: open a new issue

### Helpful links
- [Reference of all git commands](http://git-scm.com/docs).
- When you have time, I recommend reading through [the git book](http://git-scm.com/book).
- [Here's a simple guide for using git](http://rogerdudler.github.io/git-guide/), found by Suneet.
- For more details on the workflow, see
[a simple git branching model](https://gist.github.com/jbenet/ee6c9ac48068889b0912)
or
[a detailed description of the feature-branch](https://www.atlassian.com/git/workflows#!workflow-feature-branch).
