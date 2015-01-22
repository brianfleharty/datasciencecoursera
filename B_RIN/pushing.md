## How to push file to GitHub repository in windows 

###change directory to folder where your file is
	cd

###initialize a git repository in your folder
	git init

###creates connection to repository
	git remote add origin https://github.com/username/datasciencecoursera.git

###pull data from repository
	git pull origin master

###check if file exists
	ls

###list out files
	git status

###add files to staging area
	git add .

###list out files
	git status

###commits all files
	git commit -m "first commit"

###push files
	git push -u origin master


###reference guide
	http://guides.railsgirls.com/github/