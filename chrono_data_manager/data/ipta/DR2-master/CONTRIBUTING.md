## Some steps for basic contributing to IPTA DR2

1.  Make sure that you have public ssh-keys uploaded to GitLab and
have those keys available in your local ssh-agent (via ssh-add or
similar commands)

2.  [Make sure that you have your local git configured with name and email](https://git-scm.com/book/en/v2/Getting-Started-First-Time-Git-Setup)

3.  On GitLab, "Fork" the project so that you have your own copy

4.  Got to that project page and select "ssh" from the pull-down tab
so that you can clone the project locally.

5.  On your local machine, do something like the following
(substituting for YOURNAME):
```
# Clone your fork locally
> git clone git@gitlab.com:YOURNAME/DR2.git
> cd DR2
```

6.  Set up IPTA/DR2 as an "upstream" remote:
```
> git remote add upstream git@gitlab.com:IPTA/DR2.git
```

7.  Now work in your local cloned repository as needed.  You can push
your results whenever you want (using `git push`)

8.  Before you request a pull request to "upstream", make sure that
you are fully up-to-date:
```
# Make sure you are on your local master branch
> git checkout master
# Download the latest from IPTA/DR2
> git fetch upstream
# Merge those changes into your local master
> git merge upstream/master
# Alternatively, those last two commands can be combined with
> git pull upstream master
# Push everything up to your fork
> git push
```

9.  Now feel free to request a pull request to have your commits added
to IPTA/DR2.
