# install.packages("R.oo")
#library(R.oo,lib.loc="/nas/longleaf/home/weiyq/R/x86_64-redhat-linux-gnu-library/3.6")
library(R.oo)

setConstructorS3("SplitNode", function(left = -1, right = -1,curent=0, parent = +1, LP = NULL, ind = NULL, ...)
{
	extend(Object(), "SplitNode", left= left, right= right,curent=curent, parent=parent, LP=LP, ind=ind)
})

setMethodS3("ViewNode", "SplitNode", function(this,...) 
{
   return(list(left=this$left, right=this$right, parent=this$parent, LP=this$LP,ind=this$ind));
});

setMethodS3("ViewNode", "NULL", function(this,...) 
{
   return(NULL);
});

setMethodS3("get.ind", "SplitNode", function(this,...) 
{
   return(this$ind);
});


setMethodS3("get.ind", "NA", function(this,...) 
{
   return(NULL);
});
setMethodS3("get.ind", "NULL", function(this,...) 
{
   return(NULL);
});

setMethodS3("is.Node", "SplitNode", function(this,...) 
{
   return(TRUE);
});

setMethodS3("is.Node", "NULL", function(this,...) 
{
   return(FALSE);
});

setMethodS3("is.Node", "logical", function(this,...) 
{
   return(FALSE);
});

setMethodS3("noLeft", "SplitNode", function(this,...) 
{
	this$left<-NA;
});

setMethodS3("noRight", "SplitNode", function(this,...) 
{
	this$right<-NA;
});


setMethodS3("updateLP", "SplitNode", function(this,LP,...) 
{
	this$LP<-LP;
});



setMethodS3("getGroups", "SplitNode", function(this,...) 
{
	if (is.na(this$left)  && is.na(this$right) )
	{
    	return(this$curent);
  	} 
	else 
	{
    	return(c(getGroups(SplitList[[this$left]]), 
    		getGroups(SplitList[[this$right]])));  
  	}
});




groupIndex<-function(...)
{
	rtv=NULL
	grps=getGroups(SplitList[[1]])
	j=1
	for(i in grps)
	{
		rtv[j]<-list(get.ind(SplitList[[i]]))
		j=j+1
	}
	rtv
}
#grp=groupIndex()
#
#grp

#save(grp, file="test.Rdata")
setMethodS3("testSplit", "SplitNode", function(this,...) 
{
	if(!is.na(this$left) & !is.na(this$right)){ return(1) }
	else {	return(0)}
});

setMethodS3("testSplit", "NULL", function(this,...) 
{
	return(0)
});

setMethodS3("testSplit", "logical", function(this,...) 
{
	return(0)
});

stopSplit<-function(i)
{
	count=0
	for(j in (floor(i/2)+1):i)
	{
		count=count+testSplit(SplitList[[j]])
	}
	count
}


