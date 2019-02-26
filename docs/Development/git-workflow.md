# Development Workflow
The basic workflow described here follows a workflow originally [outlined](http://nvie.com/posts/a-successful-git-branching-model/) by Vincent Driessen. The workflow is built around the Git version control system. A basic description of the branching strategy and release management used by our research group is presented here. We use a central truth repository (https://github.com/NCAR/gard) that contains our main branches.

For users interested in contributing to the development, it is recommended that you contact the main developers and start an issue on the github issues page documenting the problem you wish to solve. This provides an opportunity for others to help you avoid either duplicating effort that may not yet be represented in the develop branch, or taking the code in a direction that will make it difficult to merge at a later point.

After establishing such contact, the best way to develop the code is by first "forking" the main gard repository on github so that you have your own place to develop code. You should then clone your own gard github repository to the machine you are working on locally and perform all of your development locally. When you have changes you wish to submit back to the community, push those changes to your github repository, then it the github web interface create a "pull request" to ask that these changes be pulled back into the main gard repository.

## Main Branches
Both of the main branches are published on the Github page and are controlled only by those within the admin group. The repository is organized into several branches for different purposes. In general, any new development will want to start with a branch from the develop branch.

**1. master**– The master branch represents the official release of the code. This branch is updated by new releases from the develop branch and by hotfixes.

**2. develop** – The develop branch represents the bleeding edge of the code. Because, new releases, hotfixes, and new features update the development, we recommend that all new development begins by branching from the develop branch.
## Feature Branches
For most developers, branching off the develop branch is the best choice when developing new features. This allows for easy merging of these new features into the main source code and facilitates seamless rebasing of long running feature development with the current develop branch.  We merge completed feature branches into the develop branch for inclusion in the next release.  A developer may have many feature branches, which may or may not be intended for inclusion in the main source code. Utilizing the frictionless context switching supported by Git, the development of features using this strategy provides developers with a clean and simple method of switching between features and the main branches during the development process.

## Other Branches
Although anyone could create these branches, they are designed for the preparation or hotfix of an updated master branch and therefore should only be used by members of the admin group.

**1. release** – The release branch supports the preparation of a new release. It includes any changes needed for public release or minor bug fixes.

The general workflow for creating a release is as follows:
1. Update docs/version number/etc and commit changes to master
```
# make changes
git add ...files that you changed...
git commit -m '...update version...'
```
1. Create a tag
```
git tag -a x.x.x -m 'Version x.x.x'
```
1. Push tags to github
```
git push dask master --tags
```

**2. hotfix** -  The hotfix branch facilitates mid-release bug fixes of the master branch. The key point of the hotfix branch is that it does not incorporate any new features from the develop branch, rather it is a branch off the master that addresses a specific issue or set of issues. When the hotfix is applied, the development branch is updated to reflect the hotfix changes.

## Naming Conventions
* Master branch – `master`
* Develop branch – `develop`
* Feature branch – `feature/{feature_name}`
* Hotfix branch – `hotfix/{hotfix_name}`
* Release branch – `release/{release_name}`
