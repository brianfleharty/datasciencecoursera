## How to push file to GitHub repository in windows 

#change directory to folder where your file is
	cd

#initialize a git repository in your folder
	git init

#check if file exists
	ls

#list out files
	git status

#add files to staging area
	git add .

#commits all files
	git commit -m "first commit"

#creates connection to repository
	git remote add origin https://github.com/brianfleharty/datasciencecoursera.git

#push files
	git push -u origin master



#If you want to continue making changes and pushing them to GitHub you’ll just need to 
#use the following three commands:

git add .

git commit -m "type your commit message here"

git push origin master



#reference guide
	http://guides.railsgirls.com/github/