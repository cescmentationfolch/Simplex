# Useful Git info

## Setting up a repository (already created on Github)
- Initializing repository: `git init`.
- Linking with Github: `git remote add origin https://github.com/<user>/<repo>`.
- First pull from a certain branch (e.g. `master`): `git pull origin <branch>`.

## Committing changes
- Make changes wherever and whenever you want.
- If you create a file and you want to be tracked by git, add it: `git add <file>`. If you're lazy, then `git add -A` (all files).
- Commit your change: `git commit -m <message>` or `git commit` (probably a `vi` window will appear). You will always have to write a message, so if you're lazy just write 'A'. That's not a good idea, though.
- Finally, sync with Github: `git push`. The first time, you will have to choose where to push: `git push --set-upstream origin <branch>`

## Changing between branches
- Just type `git checkout -b <branch>.

## It's not always that good
- Probably changes will be performed by others, so you will have to merge them before pushing. Before pushing, check if there are changes by `git pull`; git will automatically try to merge your changes with others'. If there is a problem, local files will be modified and you will have to decide what to do, change the file and commit and push again.
- Your last option to restore the repo from Github and erase your local files is `git fetch --all` and `git reset --hard origin/<branch>`.
