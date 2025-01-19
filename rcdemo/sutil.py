# import scipy
# import matplotlib.cm as cm
# from matplotlib.colors import LogNorm
# import matplotlib.pyplot  as plt
# import numpy as np
# import putil
# import scipy.interpolate
# import ebf
# import scipy.spatial.qhull as qhull
# import astropy.convolution

exec(open("./rcdemo/packages_to_import.py").read()) 

def interp_weights(xyz, uvw):
    d=xyz.shape[1]
    tri = qhull.Delaunay(xyz)
    simplex = tri.find_simplex(uvw)
    vertices = np.take(tri.simplices, simplex, axis=0)
    temp = np.take(tri.transform, simplex, axis=0)
    delta = uvw - temp[:, d]
    bary = np.einsum('njk,nk->nj', temp[:, :d, :], delta)
    return vertices, np.hstack((bary, 1 - bary.sum(axis=1, keepdims=True)))

def interpolate(values, vtx, wts, fill_value=np.nan):
    ret = np.einsum('nj,nj->n', np.take(values, vtx), wts)
    if fill_value is not None:
        ret[np.any(wts < 0, axis=1)] = fill_value
    return ret


def npstruct_add_field(data1,fieldname,value):
    if data1.dtype.char == 'V':
        dt=[(name,data1.dtype[name]) for name in data1.dtype.names]
        if fieldname not in data1.dtype.names:
            dt.append((fieldname,value.dtype))
            dt=np.dtype(dt)
            data2=np.zeros(data1.size,dtype=dt)
            for key in data1.dtype.names:
                data2[key]=data1[key]
            data2[fieldname]=value
            return data2
        else:
            raise RuntimeError('Field already exists:'+fieldname)
            return None
    else:
        raise RuntimeError('Data is not a numpy structure')

def chi2(x1,x2,bins=10,range=None):
    h1,be1=np.histogram(x1,range=range,bins=bins,density=False)
    h2,be2=np.histogram(x2,range=range,bins=bins,density=False)
    temp=np.sum(h1,dtype=np.float64)/np.sum(h2,dtype=np.float64)
    err=h1+h2*np.square(temp)
    ind=np.where((h1+h2)>0)[0]
    return np.sum(np.square(h1[ind]-h2[ind]*temp)/err[ind],dtype=np.float64)#/(ind.size-1)

def delta(x):
    return np.r_[x[1]-x[0],(x[2:]-x[0:-2])*0.5,x[-1]-x[-2]]

def loghist(x,bins=10,range=None,ax=None,**kwargs):
    if ax is None:
        ax1=plt
    else:
        ax1=ax
    if range is not None:
        ax1.hist(x,bins=np.logspace(np.log10(np.min(range[0])),np.log10(np.max(range[1])),bins),**kwargs)        
    else:
        ax1.hist(x,bins=np.logspace(np.log10(np.min(x)),np.log10(np.max(x)),bins),**kwargs)
    if ax is None:
        plt.xscale('log')
    else:
        ax.set_xscale('log')

def ks(x1,x2,shift=None,tol=None):
    ind=np.where(np.isfinite(x1))[0]
    x1=x1[ind]
    x1=np.sort(x1)
    y1=np.arange(x1.size)*1.0/x1.size
    if shift is not None:
        xmin=np.min(x2)*np.min(shift)*0.99
        xmax=np.max(x2)*np.max(shift)*1.01
    else:
        xmin=np.min(x2)*0.99
        xmax=np.max(x2)*1.01
    if xmin < np.min(x1):
        x1=np.r_[xmin,x1]
        y1=np.r_[0.0,y1]
    if xmax> np.max(x1):
        x1=np.r_[x1,xmax]
        y1=np.r_[y1,1.0]        
    func1=scipy.interpolate.interp1d(x1,y1)
    ind=np.where(np.isfinite(x2))[0]
    x2=x2[ind]
    x2=np.sort(x2)    
    y2=np.arange(x2.size)*1.0/x2.size
    if shift is None:
        return np.max(np.abs(func1(x2)-y2))*np.sqrt((x1.size*x2.size*1.0)/(x1.size+x2.size))
    else:
        gr=1.61803
        r1=shift[0]
        r4=shift[1]
        r2=r4-(r4-r1)/gr
        r3=r1+(r4-r1)/gr
        s1=np.max(np.abs(func1(x2*r1)-y2))
        s2=np.max(np.abs(func1(x2*r2)-y2))
        s3=np.max(np.abs(func1(x2*r3)-y2))
        s4=np.max(np.abs(func1(x2*r4)-y2))
        while np.abs(r3-r2) > tol:
            if s3 > s2:
                r4=r3
                s4=s3
                r3=r2
                s3=s2
                r2=r4-(r4-r1)/gr
                s2=np.max(np.abs(func1(x2*r2)-y2))
            else:
                r1=r2
                s1=s2
                r2=r3
                s2=s3
                r3=r1+(r4-r1)/gr
                s3=np.max(np.abs(func1(x2*r3)-y2))                
        # res=[]
        # for temp in shift:        
        #     res.append(np.max(np.abs(func1(x2*temp)-y2)))
        return r2,s2
    
#def hist2d(x, y, bins=10, range=None, normed=False, weights=None, vmin=None, vmax=None, norm=None,cmapname='jet',**kwargs):
def hist2d(x, y, bins=10, range=None, normed=False, weights=None,cmapname='jet',levels=None,clevels=None,**kwargs):
    ind=np.where(np.isfinite(x)&np.isfinite(y))[0]
    h1=hist_nd([x[ind],y[ind]],range=range,normed=normed,weights=weights,bins=bins)
    #h1.imshow(cmapname=cmapname,vmin=vmin,vmax=vmax,norm=norm,**kwargs)
    #image=h.apply(x,np.mean).reshape(h.bins)

    if levels is None:
        h1.imshow(cmapname=cmapname,dnorm=False,clevels=clevels,**kwargs)
    else:
        h1.contourf(cmapname=cmapname,levels=levels,clevels=clevels,dnorm=False,**kwargs)

    # kernel = astropy.convolution.Gaussian2DKernel(stddev=1)
    # image=h.apply(x,np.mean).reshape(h.bins)
    # temp=astropy.convolution.convolve(image.reshape(h.bins),kernel).reshape(np.prod(h.bins))
    # temp[h.data<20]=np.nan
    # h.contourf(temp,cmapname='jet', interpolation='bicubic',levels=np.linspace(-20,20,11))

def map2d(x, y, z, bins=10, range=None, normed=False, weights=None,**kwargs):
    ind=np.where(np.isfinite(x)&np.isfinite(y)&np.isfinite(z))[0]
    ind=np.where(np.isfinite(x)&np.isfinite(y)&np.isfinite(z))[0]
    h1=hist_nd([x[ind],y[ind]],range=range,normed=normed,weights=weights,bins=bins)
    #h1.imshow(cmapname=cmapname,vmin=vmin,vmax=vmax,norm=norm,**kwargs)
    image=h1.apply(z[ind],np.mean).reshape(h1.bins)
    
    h1.imshow(data=image,**kwargs)



def profile(x,z,bins=100,range=None,label=None,mincount=3,fmt='-o',color='',meanerror=False,style=None,func=None,lw=None,return_profile=False):
    """
    style: ['lines','fill_between,'errorbar']
    """
    #style=['lines','fill_between','errorbar']
    ind=np.where(np.isfinite(x)&np.isfinite(z))[0]
    x=x[ind]
    z=z[ind]
    h=hist_nd(x,bins=bins,range=range)
    ind1=np.where(h.data>=mincount)[0]
    if func is None:
        ym=h.apply(z,np.median)[ind1]
    else:
        ym=h.apply(z,func)[ind1]
    xm=h.locations[0][ind1]
    yl=h.apply(z,np.percentile,16)[ind1]
    yh=h.apply(z,np.percentile,84)[ind1]
    if meanerror:
        yl=ym+(yl-ym)/np.sqrt(h.data[ind1])
        yh=ym+(yh-ym)/np.sqrt(h.data[ind1])
    if style == 'fill_between':
        plt.fill_between(xm,yl,yh,interpolate=True,alpha=0.75)
    elif style == 'errorbar':
        plt.errorbar(xm,ym,yerr=[ym-yl,yh-ym],fmt=color+'-o',label=label,lw=lw) 
    elif style == 'lines':
        if color == '':
            color=plt.gca()._get_lines.get_next_color()
        plt.plot(xm,ym,fmt,label=label,color=color,lw=lw)
        plt.plot(xm,yl,'--',color=color,lw=lw)
        plt.plot(xm,yh,'--',color=color,lw=lw)
    else:
        plt.plot(xm,ym,color+fmt,label=label,lw=lw)
    if return_profile:
        return xm,ym

    # if errorbar:
    #     yerr=[ym+(ym-yl)/np.sqrt(h.data[ind1]),ym+(yh-ym)/np.sqrt(h.data[ind1])]
    #     plt.errorbar(xm,ym,yerr=yerr,fmt=color+'-o',label=label) 
    # elif percentile=True:
    # else:
    #     plt.plot(xm,ym,color+fmt,label=label)
    # if low is not None:
    #     yl=h.apply(z[ind],np.percentile,low)[ind1]
    #     plt.plot(xm,yl,color+'--')
    # if high is not None:
    #     yh=h.apply(z[ind],np.percentile,high)[ind1]
    #     plt.plot(xm,yh,color+'--')




class hist_nd(object):
	
	"""
	PURPOSE
	----------
		a. Makes n-dimensional histogram
		b. Returns ids for points in each bins
		c. Apply function to points in a bin
		e. Make plt.contourf and plt.imshow plots
		d. Input-> ndarray, list of ndarry, or dict of ndarray
	Copyright (c) 2015 Sanjib Sharma 
	Example:
	>>> import matplotlib.pyplot  as plt
	>>> import sutil
	>>> x=np.random.normal(size=(1000,2))
	>>> x=[np.random.normal(size=1000),np.random.normal(size=1000)]
	>>> x={'x0':np.random.normal(size=1000),'x1':np.random.normal(size=1000)}
	>>> h=sutil.hist_nd(x,bins=[10,10])
	>>> h.contourf()
	>>> h.imshow()
	>>> d=np.random.normal(size=1000)
	>>> d_mean=h.apply(d,np.mean)
	>>> h.imshow(d_mean)
	"""	
	
	def __init__(self,V,bins=10,range=None,normed=False,weights=None,keys=None,ri=True):

		"""
		Args:
	
			V: sample. The data to be histogrammed. It cab be an (N,D) array 
			   or data that can be converted to such an array, e.g.,  a list, 
			   a tuple or a dict with keys. 
	
			bins: a sequence or scalar int
	
			range: A sequence of lower and upper bin edges to be used
	
			normed: bool
	
			weights: on data points 
	
			keys: if V is disct the keys which deifne the dimensions
			  ri: bool, to specify if reverse indices are to be computed
		
		Stores:
	
			self.data: the histogram
	
			self.bins: the shape of the histogram, or no of bins along 
			each dimensions
	
			self.range: range along each dimension
	
		Methods:
			imshow: takes imshow kwargs
			colorbar: for colorbar, takes colorbar kwargs
			apply: apply function on values of 3rd variable in each bin in histdd
			resample: resample a given data to match the histogram
		"""
		self.data,self.bins,self.range,self.ri=self.compute(V,bins=bins,range=range,normed=normed,weights=weights,keys=keys,ri=ri)
	
		loc=[]
		for i in np.arange(self.bins.size):
			xe = np.linspace(self.range[i,0],self.range[i,1],self.bins[i]+1)
			xc = (xe[0:-1]+xe[1:])/2.0
			loc.append(xc)
		self.locations=loc
		self.keys=keys
	
	

	
	
	def indices(self,i=None):
		if i is None:
			return self.ri[self.ri[0]:].copy()
		if i<self.ri[0]:
			return self.ri[self.ri[i]:self.ri[i+1]]
		else:
			return np.array([],dtype='int64')



	#def contourf(self,data=None,cmapname='rainbow',levels=8,vmin=None,vmax=None,**kwargs):


	def contourf(self,data=None,cmapname='jet',levels=8,clevels=None,smooth=False,dnorm=True,**kwargs):
		"""
		Make a contour map either with input data or precomputed histogram
		Example:
		>>> h.contourf()
		or 
		>>> h.contourf(d_estimate)
		"""
		if data is None:
			data=self.data
			if dnorm:
				bs=np.float64((self.range[:,1]-self.range[:,0]))/self.bins
				ind=np.where(np.isfinite(data))[0]
				data=np.float64(data)/(np.sum(data[ind])*np.prod(bs))
				ind=np.where(data==0.0)[0]
				data[ind]=np.nan
			else:
				print( 'here')
				data=np.float64(data)
		else:
			data=data.copy()
			# ind=np.where(np.isfinite(data)==True)[0]
			# if vmin is None:
			#     vmin=np.min(data[ind])
			#     vmax=np.max(data[ind])
			#     vmin=vmin-(vmax-vmin)/255.0
			# else:
			#     ind1=ind[np.where(data[ind]<=vmin)[0]]
			#     ind2=np.where(np.isfinite(data)==False)[0]
			#     data[ind1]=vmin
			#     #data[ind1]=vmin+(vmax-vmin)/255.0
			#     #data[ind2]=vmin                
			#     ind1=ind[np.where(data[ind]>=vmax)[0]]
			#     data[ind1]=vmax                
	
		xe = np.linspace(self.range[0,0],self.range[0,1],self.bins[0]+1)
		ye = np.linspace(self.range[1,0],self.range[1,1],self.bins[1]+1)
		xc = (xe[0:-1]+xe[1:])/2.0
		yc = (ye[0:-1]+ye[1:])/2.0
		x=np.meshgrid(xc,yc,indexing='ij')
	
		data=data.reshape(self.bins)
		if smooth:
			kernel = astropy.convolution.Gaussian2DKernel(stddev=1.0)
			y=astropy.convolution.convolve(self.data.reshape(self.bins),kernel)
			data=astropy.convolution.convolve(data.reshape(self.bins),kernel)
			data[y<1]=np.nan
			kwargs['interpolation']='bicubic'
	
	
		if clevels is not None:
			plt.contour(x[0],x[1],data,levels[::len(levels)/clevels],linewidths=0.5)
			#plt.contour(x[0],x[1],data,clevels)
	
		return plt.contourf(x[0],x[1],data,levels,cmap=cm.get_cmap(cmapname),**kwargs)
        
	def contour(self,data=None,cmapname='jet',levels=8,clevels=None,smooth=False,dnorm=True,**kwargs):		
		'''
		Make a contour map either with input data or precomputed histogram
		Example:
		>>> h.contourf()
		or 
		>>> h.contourf(d_estimate)
		'''
		
		if data is None:
			data=self.data
			if dnorm:
				bs=np.float64((self.range[:,1]-self.range[:,0]))/self.bins
				ind=np.where(np.isfinite(data))[0]
				data=np.float64(data)/(np.sum(data[ind])*np.prod(bs))
				ind=np.where(data==0.0)[0]
				data[ind]=np.nan
			else:
				print( 'here')
				data=np.float64(data)
		else:
			data=data.copy()
			# ind=np.where(np.isfinite(data)==True)[0]
			# if vmin is None:
			#     vmin=np.min(data[ind])
			#     vmax=np.max(data[ind])
			#     vmin=vmin-(vmax-vmin)/255.0
			# else:
			#     ind1=ind[np.where(data[ind]<=vmin)[0]]
			#     ind2=np.where(np.isfinite(data)==False)[0]
			#     data[ind1]=vmin
			#     #data[ind1]=vmin+(vmax-vmin)/255.0
			#     #data[ind2]=vmin                
			#     ind1=ind[np.where(data[ind]>=vmax)[0]]
			#     data[ind1]=vmax                
	
		xe = np.linspace(self.range[0,0],self.range[0,1],self.bins[0]+1)
		ye = np.linspace(self.range[1,0],self.range[1,1],self.bins[1]+1)
		xc = (xe[0:-1]+xe[1:])/2.0
		yc = (ye[0:-1]+ye[1:])/2.0
		x=np.meshgrid(xc,yc,indexing='ij')
	
		data=data.reshape(self.bins)
		if smooth:
			kernel = astropy.convolution.Gaussian2DKernel(stddev=1.0)
			y=astropy.convolution.convolve(self.data.reshape(self.bins),kernel)
			data=astropy.convolution.convolve(data.reshape(self.bins),kernel)
			data[y<1]=np.nan
			kwargs['interpolation']='bicubic'
	
	
		if clevels is not None:
			plt.contour(x[0],x[1],data,levels[::len(levels)/clevels],linewidths=0.5)
			#plt.contour(x[0],x[1],data,clevels)
	
		colmap = cm.get_cmap(cmapname)
		if 'colors' in kwargs.keys():
			colmap = None
	
		return plt.contour(x[0],x[1],data,levels,cmap=colmap,**kwargs)




	def imshow(self,data=None,ax=None,cmapname='rainbow', dnorm=True,smooth=False,clevels=None,**kwargs):
		"""
		Make am image either with input data or precomputed histogram
		Example:
		>>> h.imshow()
		or 
		>>> h.imshow(d_estimate)
		"""
	#        putil.make_cmap('newjet')
		if data is None:
			data=self.data
			if dnorm:
				bs=np.float64((self.range[:,1]-self.range[:,0]))/self.bins
				ind=np.where(np.isfinite(data))[0]
				data=np.float64(data)/(np.sum(data[ind])*np.prod(bs))
				ind=np.where(data==0.0)[0]
				data[ind]=np.nan
			else:
				data=np.float64(data)
		else:
			data=data.copy()
			# ind=np.where(np.isfinite(data)==True)[0]
			# if vmin is None:
			#     vmin=np.min(data[ind])
			#     vmax=np.max(data[ind])
			#     vmin=vmin-(vmax-vmin)/255.0
			# else:
			#     ind1=ind[np.where(data[ind]<=vmin)[0]]
			#     ind2=np.where(np.isfinite(data)==False)[0]
			#     data[ind1]=vmin+(vmax-vmin)/255.0
			#     data[ind2]=vmin                
	
		d=data.reshape(self.bins).transpose()
		cmap=cm.get_cmap(cmapname)
		#cmap.set_bad('w',1.)
		#plt.imshow(np.ma.masked_where(np.isnan(d),d), interpolation=interpolation,extent=[self.range[0][0], self.range[0][1], self.range[1][0], self.range[1][1]],cmap=cmap,norm=norm,origin='lower',aspect='auto',vmin=vmin,vmax=vmax,**kwargs)
		if smooth: 
			kernel = astropy.convolution.Gaussian2DKernel(x_stddev=1.0)
			d=astropy.convolution.convolve(d,kernel)
			d[d<1]=np.nan
			kwargs['interpolation']='bicubic'
	
		xe = np.linspace(self.range[0,0],self.range[0,1],self.bins[0]+1)
		ye = np.linspace(self.range[1,0],self.range[1,1],self.bins[1]+1)
		xc = (xe[0:-1]+xe[1:])/2.0
		yc = (ye[0:-1]+ye[1:])/2.0
		x=np.meshgrid(xc,yc,indexing='ij')
	
		if ax is None:
			ax=plt
		if clevels is not None:
			ax.contour(x[0],x[1],d.transpose(),clevels,vmin=vmin,vmax=vmax)
	
		return ax.imshow(d,extent=[self.range[0][0], self.range[0][1], self.range[1][0], self.range[1][1]],cmap=cm.get_cmap(cmapname),origin='lower',aspect='auto',**kwargs)
	

	def apply(self,quant,func,*args,**kwargs):
		"""
		Apply a scalar function to data points in a bin of the histogram.
		Args:
	
			quant: an array of size N (no of data points)
	
			func: a scalar function to apply to data points          
	
		Example:
		>>> d_estimate=h.apply(d,np.mean)
			
		"""
		result=np.zeros(self.data.size,dtype=np.float64)+np.nan
		for i in np.arange(self.data.size):
			if self.data[i] > 0:
				if 'weights' in kwargs:
					result[i]=func(quant[self.indices(i)],*args,weights=kwargs['weights'][self.indices(i)])
				else:
					result[i]=func(quant[self.indices(i)],*args,**kwargs)
		return result
                
	def equisample(self,nsize):
		indt=[]
		for i in np.arange(self.data.size):
			ind=self.indices(i)
			np.random.shuffle(ind)
			if ind.size>nsize:
				indt.append(ind[0:nsize])
			else:
				indt.append(ind)
		return np.concatenate(indt)

	def resample(self,X,fsample=1.0,keys=None,verbose=False):
		if type(X) == hist_nd:
			hist1=X
		else:
			hist1=hist_nd(X,bins=self.bins,range=self.range,keys=keys,ri=True)
		hist2=self
		count1=np.sum(hist1.data)
		count2=np.sum(hist2.data)
		indices=np.zeros(np.int64(count2*fsample),dtype=np.int64)
		nprev=0
		misscount=0
		for i in np.arange(hist1.data.size):
			if hist2.data[i] > 0:
				if hist1.data[i] < hist2.data[i]*fsample: 
					misscount=misscount+hist2.data[i]*fsample-hist1.data[i]
				if hist1.data[i] > 0:
					ind1=hist1.indices(i)
					np.random.shuffle(ind1)
					nsize=np.int64(fsample*hist2.data[i])
					ind=np.arange(nsize)
					indices[nprev:nsize+nprev]=ind1[ind%ind1.size]
					nprev=nsize+nprev
	
		indices=indices[0:nprev]
		if verbose:
			print( 'fsmaple=',fsample,' primary size=',np.sum(hist2.data),' secondary size=', np.sum(hist1.data))
			print( 'Points lacking in primary=',misscount)
			print( 'Total points missed      =',int(np.sum(hist2.data)*fsample)-indices.size)
		return indices 
    

	def compute(self,X,bins=10,range=None,normed=False,weights=None,keys=None,ri=True):
		
		if type(keys) == str:
			keys=[keys]
	
	
		if type(X) == list:
			dims=len(X)
			V=X            
		elif type(X) == tuple:
			dims=len(X)
			V=X            
		elif type(X) == dict:
			dims=len(keys)
			V=[X[key] for key in keys]
		elif type(X) == np.ndarray:
			if X.dtype.char == 'V':
				dims=len(keys)
				V=[X[key] for key in keys]
			else:
				if X.ndim == 2:
					dims=X.shape[1]
					V=[X[:,i] for i in np.arange(dims)]
				elif X.ndim == 1:
					dims=1
					V=[X]
				else:
					raise RuntimeError('Input data must be of ndim=1 or 2')
		
		for i in np.arange(dims):
			if type(V[i]) != np.ndarray:
				V[i]=np.array(V[i])
	
		nsize=V[0].size
		for i in np.arange(dims):
			if V[i].size != nsize:
				print( 'dim=',i,' size=',V[i].size,' required size=',nsize)
				raise RuntimeError('all dimensions not of same size')
		#print 'dims=',dims
	
		# V={}
		# if keys==None:
		#     if X.ndim == 2:
		#         dims=X.shape[1]
		#         for i in np.arange(dims):
		#             V[i]=X[:,i]
		#     elif X.ndim == 1:
		#         dims=1
		#         V[0]=X
		#     else:
		#         raise RuntimeError('Input data must be of ndim=1 or 2')
		# else:
		#     if type(keys) == str:
		#         keys=[keys]
		#     dims=len(keys)
		#     for i,key in enumerate(keys):
		#         V[i]=X[key]
	
		#determine range and bin size
		flag=np.zeros(V[0].size,dtype=np.int64)
		nbins=np.zeros(dims,dtype=np.int64)
		nbins[:]=bins
		total_bins=np.prod(nbins)
		if range is None:
			range=np.zeros((dims,2),dtype=np.float64)
			for i in np.arange(dims):
				range[i,0]=np.min(V[i])
				range[i,1]=np.max(V[i])                
		else:
			range=np.float64(range)
			if len(range.shape) == 1:
				range=range.reshape((1,2))                    
			# clip the data to fit in range
			for i in np.arange(dims):
				# To prevent last bin from taking extra points, on 22 April 2017 changed to from 
				# ind=np.where((V[i]<range[i,0])|(V[i]>range[i,1]))[0]
				ind=np.where((V[i]<range[i,0])|(V[i]>=range[i,1]))[0]
				flag[ind]=1                
		bs=np.float64((range[:,1]-range[:,0]))/nbins
		ind=np.where(flag==0)[0]
	
		#index of each data point
		d=np.clip(np.int64((V[0]-range[0,0])/bs[0]),0,nbins[0]-1)
		for i in np.arange(1,dims):                     
			# The scaled indices, 
			d=nbins[i]*d+np.clip(np.int64((V[i]-range[i,0])/bs[i]),0,nbins[i]-1)
		d=d[ind]
	
	#        h,e=np.histogram(d,range=[0,total_bins],bins=total_bins,normed=normed,weights=weights)
	
		# np.bincount is faster
		h=np.bincount(d,weights=None,minlength=total_bins)
		e=np.arange(total_bins+1)
	
		# reverse indices
		if ri == True:
			inds=np.argsort(d)
			rev_ind=np.zeros(e.size+inds.size,dtype='int64') 
			rev_ind[0]=e.size        
			rev_ind[1:e.size]=np.cumsum(h)+e.size        
			rev_ind[e.size:]=ind[inds]
		else:
			rev_ind=None
	
		if weights is not None:
			h=np.bincount(d,weights=weights[ind],minlength=total_bins)
	
		if normed==True:
			h=np.float64(h)/(np.sum(h)*np.prod(bs))
	
		return h,nbins,range,rev_ind
	



class hist_ndp(hist_nd):    
    def imshow1(self,density=False,vmin=None,vmax=None,**kwargs):
#        putil.make_cmap('newjet')
        data=self.data.copy()
        if density:
            bs=np.float64((self.range[:,1]-self.range[:,0]))/self.bins
            data=np.float64(data)/(np.sum(data)*np.prod(bs))
        else:
            ind=np.where(np.isfinite(data)==True)[0]
            if vmin is None:
                vmin=np.min(data[ind])
                vmax=np.max(data[ind])
                vmin=vmin-(vmax-vmin)/255.0
            else:
                ind1=ind[np.where(data[ind]<=vmin)[0]]
                ind2=np.where(np.isfinite(data)==False)[0]
                data[ind1]=vmin+(vmax-vmin)/255.0
                data[ind2]=vmin                
        plt.imshow(data.reshape(self.bins).transpose(),extent=[self.range[0][0], self.range[0][1], self.range[1][0], self.range[1][1]],vmin=vmin,vmax=vmax,origin='lower',aspect='auto',**kwargs)
    # def recbin(x, y, C=None, gridsize=100, bins=None, xscale='linear', yscale='linear', extent=None, cmap=None, norm=None, vmin=None, vmax=None, alpha=None, linewidths=None, edgecolors='none', reduce_C_function=<function mean>, mincnt=None, marginals=False, hold=None, **kwargs):
    # h1=hist_nd([x,y],range=range,normed=normed,weights=weights)



			
		








class StructGrid():
    """
    Takes an ndarray of npstruct and makes a grid out of it 
    for use in interpolation.
    """
    @staticmethod
    def groupby(data,keys,regular=True):
        """
        Will group data 
        Parameters
        ----------
        data :  (ndarray of npstruct)
        keys : a list of strings specifying struct field names 
               to group by
        Returns
        -------
        res   : (list) data split in groups
        groups: (list) of length len(keys)
                 groups[j][i] value of key[j] in i-th group         
        """
        temp=0
        for key in keys:
            temp=temp+np.abs(np.diff(data[key]))
        res=np.split(data,np.where(temp!=0)[0]+1)

        groups1=[]
        for key in keys:
            groups1.append([])
        for temp in res:
            for i,key in enumerate(keys):
                groups1[i].append(temp[key][0])
        groups1=[np.array(group) for group in groups1]

        if regular==False:
            return res,groups1
        else:
            groups2=[np.unique(data[key]) for key in keys]
            myshape=[group.size for group in groups2]
            points=(np.meshgrid(*groups2,indexing='ij'))
            groups2=[np.ravel(point) for point in points]


            print( 'groupby keys:',keys)
            print( 'group size:',groups1[0].size)
            print( 'grid size :',groups2[0].size,', ',myshape )

            j=0
            for i in range(groups2[0].size):
                temp=0.0
                for k in range(len(groups1)):
                    temp=temp+np.abs(groups1[k][j]-groups2[k][i])
                if temp != 0:
                    res.insert(i,None)
                else:
                    j=j+1

            return [res,groups2]

    @staticmethod
    def fill_missing(data,keys,miss_order=None):        
        """
        
        Parameters
        ----------
        data :  (ndarray of npstruct)
        keys : a list of strings specifying struct field names 
               to group by
        Returns
        -------
        res   : (list) data split in groups
        groups: (list) of length len(keys)
                 groups[j][i] value of key[j] in i-th group         
        """
        data.sort(order=keys)    
        vnames=[temp for temp in data.dtype.names if temp not in keys]
        myshape=np.array([np.unique(data[key]).size for key in keys])
        ind=np.where(np.isfinite(data[vnames[0]])==False)[0]
        if miss_order is None:
            miss_order=range(len(keys))
        else:
            if type(miss_order[0]) == str:
                miss_order=[keys.index(temp) for temp in miss_order]
        print( 'Filling missing grid points------>' )
        print( 'grid shape:',myshape)
        print( 'missing ',ind.size,' out of ', data.size )
        idone=0
        idone1=0
        for i in ind:
            index=np.unravel_index(i,myshape)
            ibk=None
            for j in miss_order:
                k=index[j]
                i1=i
                i2=i
                while (k+1)<myshape[j]:
                    k=k+1
                    index1=list(index)
                    index1[j]=k
                    temp=np.ravel_multi_index(index1,myshape)
                    if np.isfinite(data[vnames[0]][temp]) == True:
                        i1=temp
                        break
                k=index[j]
                while (k-1)>0:
                    k=k-1
                    index1=list(index)
                    index1[j]=k
                    temp=np.ravel_multi_index(index1,myshape)
                    if np.isfinite(data[vnames[0]][temp]) == True:
                        i2=temp
                        break
                if (i1!=i)and(i2!=i):
                    idone=idone+1
                    u=(data[keys[j]][i]-data[keys[j]][i1])/(data[keys[j]][i2]-data[keys[j]][i1])
                    for vname in vnames:
                        data[vname][i]=data[vname][i1]*(1-u)+data[vname][i2]*u
                    break
                if ibk is None:
                    if i1!=i:
                        ibk=i1
                    elif i2 != i:
                        ibk=i2
            if np.isfinite(data[vnames[0]][i]) == False:
                if ibk is not None:
                    for vname in vnames:
                        data[vname][i]=data[vname][ibk]
                    idone=idone+1
                    idone1=idone1+1

        print(idone  ,' out of ', ind.size,' filled (', idone1, ' using nearest)')

    @staticmethod
    def griddata(data1,keys,qname,quant,miss_order=None):
        """
        Construct a regular grid out of given points. 
        One dimension can be regridded if it is not already regular.
        If few grid points are missing it will also fill them. 
        Parameters
        ----------
        data1 : (ndarray of npstruct)
        keys  : a list of strings specifying struct field names 
                which will be grid dimensions
        qname : field name which is irregular and needs to be regridded
                qname should be present in list keys            
        quant : the array specifying new values of qname
        Returns
        -------
        data_new: (ndarray of npstruct) conforming to a regular grid,  
                  sorted and ordered by keys.
        """
        data=np.resize(data1,data1.size)
        if qname not in keys:
            raise RuntimeError('qname should be present in keys')
        split_keys=[temp for temp in keys if temp != qname]
        vnames=[temp for temp in data.dtype.names if temp not in keys]
        data.sort(order=split_keys+[qname])
        data_split,groups=StructGrid.groupby(data,keys=split_keys)

        mylist=[]
        for j,d in enumerate(data_split):
            data2=np.zeros(quant.size,dtype=data.dtype)
            data2[qname]=quant
            for i,key in enumerate(split_keys):
                data2[key]=groups[i][j]
            if d is None:
                for vname in vnames:
                    data2[vname]=np.nan
            else:
                ind1=np.argsort(d[qname])                    
                xs=np.r_[quant[0],d[qname][ind1],quant[-1]]
                for i,vname in enumerate(vnames):
                    ys=np.r_[d[vname][ind1[0]],d[vname][ind1],d[vname][ind1[-1]]]
                    data2[vname]=scipy.interpolate.interp1d(xs,ys,assume_sorted=True)(data2[qname])
            mylist.append(data2)
                

        data_new=np.concatenate(mylist)
        StructGrid.fill_missing(data_new,keys,miss_order)
        data_new.sort(order=keys)
        return data_new








def make_RegularGridInterpolator(outfile,gridname,dimnames,points,myfunc,**kwargs):
    bins=[temp.size for temp in points]
    points1=np.meshgrid(*points,indexing='ij')
    points1=np.array([temp.reshape(-1) for temp in points1]).transpose()
    temp=np.zeros(points1[:,0].size,dtype=np.float64)
    for i in range(points1.shape[0]):
        if i % 100 == 0:
            print( i,' out of ',temp.size)
        temp[i]=myfunc(*points1[i,:],**kwargs)
    values=temp.reshape(bins)

    if outfile is not None:
        if gridname[0].startswith('/') == False:
            gridname='/'+gridname
        ebf.write(outfile,gridname+'/data',values,'w')
        ebf.write(outfile,gridname+'/dimnames',dimnames,'a')
        for i,temp in enumerate(points):
            ebf.write(outfile,gridname+'/'+dimnames[i],points[i],'a')

    return scipy.interpolate.RegularGridInterpolator(points,values,bounds_error=False,fill_value=np.nan)

def read_RegularGridInterpolator(outfile,gridname):
    if gridname[0].startswith('/') == False:
        gridname='/'+gridname
    values=ebf.read(outfile,gridname+'/data')
    dimnames=ebf.read(outfile,gridname+'/dimnames')
    points=[ebf.read(outfile,gridname+'/'+dimname) for dimname in dimnames]
    return scipy.interpolate.RegularGridInterpolator(points,values,bounds_error=False,fill_value=np.nan)

class InterpGrid():
    """
    Takes an ndarray of npstruct and makes a grid out of it 
    for use in interpolation.
    """
    def __init__(self,data1,keys):
        """
        Parameters
        ----------
        data1 :  (ndarray of npstruct)
        keys : a list of strings specifying struct field names 
               to group by
        """
        data=np.resize(data1,data1.size)
        data.sort(order=keys)
        self.keys=keys
        self.vnames=[temp for temp in data.dtype.names if temp not in self.keys]

        self.points=[np.unique(data[key]) for key in self.keys]
        self.delta=[self._get_delta(x) for x in self.points]

        self.values={}; print(self.points)        
        for vname in self.vnames:
            self.values[vname]=data[vname].reshape([point.size for point in self.points])

        self.points1=tuple([data[key] for key in self.keys])
        self.values1={}        
        temp=np.meshgrid(*self.delta,indexing='ij')
        self.vol=1.0
        for x in temp:
            self.vol=self.vol*x
        for vname in self.vnames:
            self.values1[vname]=data[vname]
    
    def _get_delta(self,x):
        return np.r_[x[1]-x[0],(x[0:-2]-x[2:])*0.5,x[-1]-x[-2]]

    def _homogenize_arrays(self,xi):
        xj=[np.asarray(t) for t in xi]
        temp=xj[0]
        for t in xj:
            temp=temp+t
        xj=[np.zeros_like(temp)+t for t in xj]
        return xj

    def get_values(self,vname,xi,fill_value=None):
        """
        Parameters
        ----------
        vname: field name whose interpolated value is desired
        xi   : the coordinates for which interpolation is required
               should be in same order as keys
        Returns
        -------
        value: the interpolated value
        """
        fill_value1=np.nan
        if type(xi) == list:
            xi=np.array(self._homogenize_arrays(xi)).transpose()
        t1=scipy.interpolate.interpn(self.points,self.values[vname],xi,bounds_error=False,fill_value=fill_value1)
        if fill_value == 'nearest':
            ind=np.where(np.isfinite(t1)==False)[0]
            if ind.size>0:
                print( 'outside interp range',ind.size,' out of ',t1.size )
                t1[ind]=scipy.interpolate.griddata(self.points1,self.values1[vname],xi[ind],method='nearest') 
        return t1




##def old_sutil():
	
	###---------- scripts from Sanjib:-----------------------------------------------------------------------------------#
	
	
	##class hist_nd_old(object):
	##'''
	##PURPOSE
	##----------
		##a. Makes n-dimensional histogram
		##b. Returns ids for points in each bins
		##c. Apply function to points in a bin
		##e. Make plt.contourf and plt.imshow plots
		##d. Input-> ndarray, list of ndarry, or dict of ndarray
	##Copyright (c) 2015 Sanjib Sharma 
	##Example:
	##>>> import matplotlib.pyplot  as plt
	##>>> import sutil
	##>>> x=np.random.normal(size=(1000,2))
	##>>> x=[np.random.normal(size=1000),np.random.normal(size=1000)]
	##>>> x={'x0':np.random.normal(size=1000),'x1':np.random.normal(size=1000)}
	##>>> h=sutil.hist_nd(x,bins=[10,10])
	##>>> h.contourf()
	##>>> h.imshow()
	##>>> d=np.random.normal(size=1000)
	##>>> d_mean=h.apply(d,np.mean)
	##>>> h.imshow(d_mean)
	##'''

	##def __init__(self,V,bins=10,range=None,normed=False,weights=None,keys=None,ri=True):
		##"""
		##Args:
	
			##V: sample. The data to be histogrammed. It cab be an (N,D) array 
			   ##or data that can be converted to such an array, e.g.,  a list, 
			   ##a tuple or a dict with keys. 
	
			##bins: a sequence or scalar int
	
			##range: A sequence of lower and upper bin edges to be used
	
			##normed: bool
	
			##weights: on data points 
	
			##keys: if V is disct the keys which deifne the dimensions
			  ##ri: bool, to specify if reverse indices are to be computed
		
		##Stores:
	
			##self.data: the histogram
	
			##self.bins: the shape of the histogram, or no of bins along 
			##each dimensions
	
			##self.range: range along each dimension
	
		##Methods:
			##imshow: takes imshow kwargs
			##colorbar: for colorbar, takes colorbar kwargs
			##apply: apply function on values of 3rd variable in each bin in histdd
			##resample: resample a given data to match the histogram
		##"""
		##self.data,self.bins,self.range,self.ri=self.compute(V,bins=bins,range=range,normed=normed,weights=weights,keys=keys,ri=ri)
		##loc=[]
		##for i in np.arange(self.bins.size):
			##xe = np.linspace(self.range[i,0],self.range[i,1],self.bins[i]+1)
			##xc = (xe[0:-1]+xe[1:])/2.0
			##loc.append(xc)
		##self.locations=loc
		##self.keys=keys
	
	##def indices(self,i=None):
		##if i is None:
			##return self.ri[self.ri[0]:].copy()
		##if i<self.ri[0]:
			##return self.ri[self.ri[i]:self.ri[i+1]]
		##else:
			##return np.array([],dtype='int64')
	
	###def contourf(self,data=None,cmapname='rainbow',levels=8,vmin=None,vmax=None,**kwargs):
	##def contourf(self,data=None,cmapname='jet',levels=8,clevels=None,smooth=False,kernel_dev=1.0,dnorm=True,vmin=None,vmax=None,**kwargs):
		##"""
		##Make a contour map either with input data or precomputed histogram
		##Example:
		##>>> h.contourf()
		##or 
		##>>> h.contourf(d_estimate)
		##"""
		##if data is None:
			##data=self.data
			##if dnorm:
				##bs=np.float64((self.range[:,1]-self.range[:,0]))/self.bins
				##ind=np.where(np.isfinite(data))[0]
				##data=np.float64(data)/(np.sum(data[ind])*np.prod(bs))
				##ind=np.where(data==0.0)[0]
				##data[ind]=np.nan
			##else:
				##print 'here'
				##data=np.float64(data)
		##else:
			##data=data.copy()
			### ind=np.where(np.isfinite(data)==True)[0]
			### if vmin is None:
			###     vmin=np.min(data[ind])
			###     vmax=np.max(data[ind])
			###     vmin=vmin-(vmax-vmin)/255.0
			### else:
			###     ind1=ind[np.where(data[ind]<=vmin)[0]]
			###     ind2=np.where(np.isfinite(data)==False)[0]
			###     data[ind1]=vmin
			###     #data[ind1]=vmin+(vmax-vmin)/255.0
			###     #data[ind2]=vmin                
			###     ind1=ind[np.where(data[ind]>=vmax)[0]]
			###     data[ind1]=vmax                
	
		##xe = np.linspace(self.range[0,0],self.range[0,1],self.bins[0]+1)
		##ye = np.linspace(self.range[1,0],self.range[1,1],self.bins[1]+1)
		##xc = (xe[0:-1]+xe[1:])/2.0
		##yc = (ye[0:-1]+ye[1:])/2.0
		##x=np.meshgrid(xc,yc,indexing='ij')
	
		##data=data.reshape(self.bins)
		##if smooth:
			##from astropy.convolution import Gaussian2DKernel
			##kernel = Gaussian2DKernel(stddev=kernel_dev)
			##y=astropy.convolution.convolve(self.data.reshape(self.bins),kernel)
			##data=astropy.convolution.convolve(data.reshape(self.bins),kernel)
			##data[y<1]=np.nan
			##kwargs['interpolation']='bicubic'
	
	
		##if clevels is not None:
			##plt.contour(x[0],x[1],data,levels[::len(levels)/clevels],linewidths=0.5)
			###plt.contour(x[0],x[1],data,clevels)
	
		##plt.contourf(x[0],x[1],data,levels,cmap=cm.get_cmap(cmapname),vmin=vmin,vmax=vmax,**kwargs)
	
	
	##def imshow(self,data=None,cmapname='rainbow', dnorm=True,smooth=False,kernel_dev=1.0,clevels=None,vmin=None,vmax=None,**kwargs):
		##"""
		##Make am image either with input data or precomputed histogram
		##Example:
		##>>> h.imshow()
		##or 
		##>>> h.imshow(d_estimate)
		##"""
	###        putil.make_cmap('newjet')
		##if data is None:
			##data=self.data
			##if dnorm:
				##bs=np.float64((self.range[:,1]-self.range[:,0]))/self.bins
				##ind=np.where(np.isfinite(data))[0]
				##data=np.float64(data)/(np.sum(data[ind])*np.prod(bs))
				##ind=np.where(data==0.0)[0]
				##data[ind]=np.nan
			##else:
				##data=np.float64(data)
		##else:
			##data=data.copy()
			### ind=np.where(np.isfinite(data)==True)[0]
			### if vmin is None:
			###     vmin=np.min(data[ind])
			###     vmax=np.max(data[ind])
			###     vmin=vmin-(vmax-vmin)/255.0
			### else:
			###     ind1=ind[np.where(data[ind]<=vmin)[0]]
			###     ind2=np.where(np.isfinite(data)==False)[0]
			###     data[ind1]=vmin+(vmax-vmin)/255.0
			###     data[ind2]=vmin                
	
		##d=data.reshape(self.bins).transpose()
		##cmap=cm.get_cmap(cmapname)
		###cmap.set_bad('w',1.)
		###plt.imshow(np.ma.masked_where(np.isnan(d),d), interpolation=interpolation,extent=[self.range[0][0], self.range[0][1], self.range[1][0], self.range[1][1]],cmap=cmap,norm=norm,origin='lower',aspect='auto',vmin=vmin,vmax=vmax,**kwargs)
		##if smooth: 
			##kernel = astropy.convolution.Gaussian2DKernel(stddev=kernel_dev)
			##d=astropy.convolution.convolve(d,kernel)
			###d[d<1]=np.nan
			##kwargs['interpolation']= 'bicubic' # 'nearest'
	
		##xe = np.linspace(self.range[0,0],self.range[0,1],self.bins[0]+1)
		##ye = np.linspace(self.range[1,0],self.range[1,1],self.bins[1]+1)
		##xc = (xe[0:-1]+xe[1:])/2.0
		##yc = (ye[0:-1]+ye[1:])/2.0
		##x=np.meshgrid(xc,yc,indexing='ij')
	
		##if clevels is not None:
			##plt.contour(x[0],x[1],d.transpose(),clevels,vmin=vmin,vmax=vmax)
	
		##plt.imshow(d,extent=[self.range[0][0], self.range[0][1], self.range[1][0], self.range[1][1]],cmap=cm.get_cmap(cmapname),origin='lower',aspect='auto',vmin=vmin, vmax=vmax,**kwargs)
	
	##def apply(self,quant,func,*args,**kwargs):
		##"""
		##Apply a scalar function to data points in a bin of the histogram.
		##Args:
	
			##quant: an array of size N (no of data points)
	
			##func: a scalar function to apply to data points          
	
		##Example:
		##>>> d_estimate=h.apply(d,np.mean)
			
		##"""
		##result=np.zeros(self.data.size,dtype=np.float64)+np.nan
		##for i in np.arange(self.data.size):
			##if self.data[i] > 0:
				##if 'weights' in kwargs:
					##result[i]=func(quant[self.indices(i)],*args,weights=kwargs['weights'][self.indices(i)])
				##else:
					##result[i]=func(quant[self.indices(i)],*args,**kwargs)
		##return result
	##def equisample(self,nsize):
		##indt=[]
		##for i in np.arange(self.data.size):
			##ind=self.indices(i)
			##np.random.shuffle(ind)
			##if ind.size>nsize:
				##indt.append(ind[0:nsize])
			##else:
				##indt.append(ind)
		##return np.concatenate(indt)
	##def resample(self,X,fsample=1.0,keys=None,verbose=False):
		##if type(X) == hist_nd:
			##hist1=X
		##else:
			##hist1=hist_nd(X,bins=self.bins,range=self.range,keys=keys,ri=True)
		##hist2=self
		##count1=np.sum(hist1.data)
		##count2=np.sum(hist2.data)
		##indices=np.zeros(np.int64(count2*fsample),dtype=np.int64)
		##nprev=0
		##misscount=0
		##for i in np.arange(hist1.data.size):
			##if hist2.data[i] > 0:
				##if hist1.data[i] < hist2.data[i]*fsample: 
					##misscount=misscount+hist2.data[i]*fsample-hist1.data[i]
				##if hist1.data[i] > 0:
					##ind1=hist1.indices(i)
					##np.random.shuffle(ind1)
					##nsize=np.int64(fsample*hist2.data[i])
					##ind=np.arange(nsize)
					##indices[nprev:nsize+nprev]=ind1[ind%ind1.size]
					##nprev=nsize+nprev
	
		##indices=indices[0:nprev]
		##if verbose:
			##print 'fsmaple=',fsample,' primary size=',np.sum(hist2.data),' secondary size=', np.sum(hist1.data)
			##print 'Points lacking in primary=',misscount
			##print 'Total points missed      =',int(np.sum(hist2.data)*fsample)-indices.size
		##return indices 
	##def compute(self,X,bins=10,range=None,normed=False,weights=None,keys=None,ri=True):
		##if type(keys) == str:
			##keys=[keys]
	
		##if type(X) == list:
			##dims=len(X)
			##V=X            
		##elif type(X) == tuple:
			##dims=len(X)
			##V=X            
		##elif type(X) == dict:
			##dims=len(keys)
			##V=[X[key] for key in keys]
		##elif type(X) == np.ndarray:
			##if X.dtype.char == 'V':
				##dims=len(keys)
				##V=[X[key] for key in keys]
			##else:
				##if X.ndim == 2:
					##dims=X.shape[1]
					##V=[X[:,i] for i in np.arange(dims)]
				##elif X.ndim == 1:
					##dims=1
					##V=[X]
				##else:
					##raise RuntimeError('Input data must be of ndim=1 or 2')
		
		##for i in np.arange(dims):
			##if type(V[i]) != np.ndarray:
				##V[i]=np.array(V[i])
	
		##nsize=V[0].size
		##for i in np.arange(dims):
			##if V[i].size != nsize:
				##print 'dim=',i,' size=',V[i].size,' required size=',nsize
				##raise RuntimeError('all dimensions not of same size')
		###print 'dims=',dims
	
		### V={}
		### if keys==None:
		###     if X.ndim == 2:
		###         dims=X.shape[1]
		###         for i in np.arange(dims):
		###             V[i]=X[:,i]
		###     elif X.ndim == 1:
		###         dims=1
		###         V[0]=X
		###     else:
		###         raise RuntimeError('Input data must be of ndim=1 or 2')
		### else:
		###     if type(keys) == str:
		###         keys=[keys]
		###     dims=len(keys)
		###     for i,key in enumerate(keys):
		###         V[i]=X[key]
	
		###determine range and bin size
		##flag=np.zeros(V[0].size,dtype=np.int64)
		##nbins=np.zeros(dims,dtype=np.int64)
		##nbins[:]=bins
		##total_bins=np.prod(nbins)
		##if range is None:
			##range=np.zeros((dims,2),dtype=np.float64)
			##for i in np.arange(dims):
				##range[i,0]=np.min(V[i])
				##range[i,1]=np.max(V[i])                
		##else:
			##range=np.float64(range)
			##if len(range.shape) == 1:
				##range=range.reshape((1,2))                    
			### clip the data to fit in range
			##for i in np.arange(dims):
				### To prevent last bin from taking extra points, on 22 April 2017 changed to from 
				### ind=np.where((V[i]<range[i,0])|(V[i]>range[i,1]))[0]
				##ind=np.where((V[i]<range[i,0])|(V[i]>=range[i,1]))[0]
				##flag[ind]=1                
		##bs=np.float64((range[:,1]-range[:,0]))/nbins
		##ind=np.where(flag==0)[0]
	
		###index of each data point
		##d=np.clip(np.int64((V[0]-range[0,0])/bs[0]),0,nbins[0]-1)
		##for i in np.arange(1,dims):                     
			### The scaled indices, 
			##d=nbins[i]*d+np.clip(np.int64((V[i]-range[i,0])/bs[i]),0,nbins[i]-1)
		##d=d[ind]
	
	###        h,e=np.histogram(d,range=[0,total_bins],bins=total_bins,normed=normed,weights=weights)
	
		### np.bincount is faster
		##h=np.bincount(d,weights=None,minlength=total_bins)
		##e=np.arange(total_bins+1)
	
		### reverse indices
		##if ri == True:
			##inds=np.argsort(d)
			##rev_ind=np.zeros(e.size+inds.size,dtype='int64') 
			##rev_ind[0]=e.size        
			##rev_ind[1:e.size]=np.cumsum(h)+e.size        
			##rev_ind[e.size:]=ind[inds]
		##else:
			##rev_ind=None
	
		##if weights is not None:
			##h=np.bincount(d,weights=weights[ind],minlength=total_bins)
	
		##if normed==True:
			##h=np.float64(h)/(np.sum(h)*np.prod(bs))
	
		##return h,nbins,range,rev_ind
	
	
	##class hist_nd(object):
	##"""
	##PURPOSE
	##----------
		##a. Makes n-dimensional histogram
		##b. Returns ids for points in each bins
		##c. Apply function to points in a bin
		##e. Make plt.contourf and plt.imshow plots
		##d. Input-> ndarray, list of ndarry, or dict of ndarray
	##Copyright (c) 2015 Sanjib Sharma 
	##Example:
	##>>> import matplotlib.pyplot  as plt
	##>>> import sutil
	##>>> x=np.random.normal(size=(1000,2))
	##>>> x=[np.random.normal(size=1000),np.random.normal(size=1000)]
	##>>> x={'x0':np.random.normal(size=1000),'x1':np.random.normal(size=1000)}
	##>>> h=sutil.hist_nd(x,bins=[10,10])
	##>>> h.contourf()
	##>>> h.imshow()
	##>>> d=np.random.normal(size=1000)
	##>>> d_mean=h.apply(d,np.mean)
	##>>> h.imshow(d_mean)
	##"""
	##def __init__(self,V,bins=10,range=None,normed=False,weights=None,keys=None,ri=True):
		##"""
		##Args:
	
			##V: sample. The data to be histogrammed. It cab be an (N,D) array 
			   ##or data that can be converted to such an array, e.g.,  a list, 
			   ##a tuple or a dict with keys. 
	
			##bins: a sequence or scalar int
	
			##range: A sequence of lower and upper bin edges to be used
	
			##normed: bool
	
			##weights: on data points 
	
			##keys: if V is disct the keys which deifne the dimensions
			  ##ri: bool, to specify if reverse indices are to be computed
		
		##Stores:
	
			##self.data: the histogram
	
			##self.bins: the shape of the histogram, or no of bins along 
			##each dimensions
	
			##self.range: range along each dimension
	
		##Methods:
			##imshow: takes imshow kwargs
			##colorbar: for colorbar, takes colorbar kwargs
			##apply: apply function on values of 3rd variable in each bin in histdd
			##resample: resample a given data to match the histogram
		##"""
		##self.data,self.bins,self.range,self.ri=self.compute(V,bins=bins,range=range,normed=normed,weights=weights,keys=keys,ri=ri)
		##loc=[]
		##for i in np.arange(self.bins.size):
			##xe = np.linspace(self.range[i,0],self.range[i,1],self.bins[i]+1)
			##xc = (xe[0:-1]+xe[1:])/2.0
			##loc.append(xc)
		##self.locations=loc
		##self.keys=keys
	##def indices(self,i=None):
		##if i is None:
			##return self.ri[self.ri[0]:].copy()
		##if i<self.ri[0]:
			##return self.ri[self.ri[i]:self.ri[i+1]]
		##else:
			##return np.array([],dtype='int64')
	###def contourf(self,data=None,cmapname='rainbow',levels=8,vmin=None,vmax=None,**kwargs):
	##def contourf(self,data=None,cmapname='jet',levels=8,clevels=None,smooth=False,dnorm=True,**kwargs):
		##"""
		##Make a contour map either with input data or precomputed histogram
		##Example:
		##>>> h.contourf()
		##or 
		##>>> h.contourf(d_estimate)
		##"""
		##if data is None:
			##data=self.data
			##if dnorm:
				##bs=np.float64((self.range[:,1]-self.range[:,0]))/self.bins
				##ind=np.where(np.isfinite(data))[0]
				##data=np.float64(data)/(np.sum(data[ind])*np.prod(bs))
				##ind=np.where(data==0.0)[0]
				##data[ind]=np.nan
			##else:
				##print 'here'
				##data=np.float64(data)
		##else:
			##data=data.copy()
			### ind=np.where(np.isfinite(data)==True)[0]
			### if vmin is None:
			###     vmin=np.min(data[ind])
			###     vmax=np.max(data[ind])
			###     vmin=vmin-(vmax-vmin)/255.0
			### else:
			###     ind1=ind[np.where(data[ind]<=vmin)[0]]
			###     ind2=np.where(np.isfinite(data)==False)[0]
			###     data[ind1]=vmin
			###     #data[ind1]=vmin+(vmax-vmin)/255.0
			###     #data[ind2]=vmin                
			###     ind1=ind[np.where(data[ind]>=vmax)[0]]
			###     data[ind1]=vmax                
	
		##xe = np.linspace(self.range[0,0],self.range[0,1],self.bins[0]+1)
		##ye = np.linspace(self.range[1,0],self.range[1,1],self.bins[1]+1)
		##xc = (xe[0:-1]+xe[1:])/2.0
		##yc = (ye[0:-1]+ye[1:])/2.0
		##x=np.meshgrid(xc,yc,indexing='ij')
	
		##data=data.reshape(self.bins)
		##if smooth:
			##kernel = astropy.convolution.Gaussian2DKernel(stddev=1.0)
			##y=astropy.convolution.convolve(self.data.reshape(self.bins),kernel)
			##data=astropy.convolution.convolve(data.reshape(self.bins),kernel)
			##data[y<1]=np.nan
			##kwargs['interpolation']='bicubic'
	
	
		##if clevels is not None:
			##plt.contour(x[0],x[1],data,levels[::len(levels)/clevels],linewidths=0.5)
			###plt.contour(x[0],x[1],data,clevels)
	
		##return plt.contourf(x[0],x[1],data,levels,cmap=cm.get_cmap(cmapname),**kwargs)
	##def contour(self,data=None,cmapname='jet',levels=8,clevels=None,smooth=False,dnorm=True,**kwargs):
		##"""
		##Make a contour map either with input data or precomputed histogram
		##Example:
		##>>> h.contourf()
		##or 
		##>>> h.contourf(d_estimate)
		##"""
		##if data is None:
			##data=self.data
			##if dnorm:
				##bs=np.float64((self.range[:,1]-self.range[:,0]))/self.bins
				##ind=np.where(np.isfinite(data))[0]
				##data=np.float64(data)/(np.sum(data[ind])*np.prod(bs))
				##ind=np.where(data==0.0)[0]
				##data[ind]=np.nan
			##else:
				##print 'here'
				##data=np.float64(data)
		##else:
			##data=data.copy()
			### ind=np.where(np.isfinite(data)==True)[0]
			### if vmin is None:
			###     vmin=np.min(data[ind])
			###     vmax=np.max(data[ind])
			###     vmin=vmin-(vmax-vmin)/255.0
			### else:
			###     ind1=ind[np.where(data[ind]<=vmin)[0]]
			###     ind2=np.where(np.isfinite(data)==False)[0]
			###     data[ind1]=vmin
			###     #data[ind1]=vmin+(vmax-vmin)/255.0
			###     #data[ind2]=vmin                
			###     ind1=ind[np.where(data[ind]>=vmax)[0]]
			###     data[ind1]=vmax                
	
		##xe = np.linspace(self.range[0,0],self.range[0,1],self.bins[0]+1)
		##ye = np.linspace(self.range[1,0],self.range[1,1],self.bins[1]+1)
		##xc = (xe[0:-1]+xe[1:])/2.0
		##yc = (ye[0:-1]+ye[1:])/2.0
		##x=np.meshgrid(xc,yc,indexing='ij')
	
		##data=data.reshape(self.bins)
		##if smooth:
			##kernel = astropy.convolution.Gaussian2DKernel(stddev=1.0)
			##y=astropy.convolution.convolve(self.data.reshape(self.bins),kernel)
			##data=astropy.convolution.convolve(data.reshape(self.bins),kernel)
			##data[y<1]=np.nan
			##kwargs['interpolation']='bicubic'
	
	
		##if clevels is not None:
			##plt.contour(x[0],x[1],data,levels[::len(levels)/clevels],linewidths=0.5)
			###plt.contour(x[0],x[1],data,clevels)
	
		##colmap = cm.get_cmap(cmapname)
		##if 'colors' in kwargs.keys():
			##colmap = None
	
		##return plt.contour(x[0],x[1],data,levels,cmap=colmap,**kwargs)
		
		
	##def imshow(self,data=None,ax=None,cmapname='rainbow', dnorm=True,smooth=False,clevels=None,**kwargs):
		##"""
		##Make am image either with input data or precomputed histogram
		##Example:
		##>>> h.imshow()
		##or 
		##>>> h.imshow(d_estimate)
		##"""
	###        putil.make_cmap('newjet')
		##if data is None:
			##data=self.data
			##if dnorm:
				##bs=np.float64((self.range[:,1]-self.range[:,0]))/self.bins
				##ind=np.where(np.isfinite(data))[0]
				##data=np.float64(data)/(np.sum(data[ind])*np.prod(bs))
				##ind=np.where(data==0.0)[0]
				##data[ind]=np.nan
			##else:
				##data=np.float64(data)
		##else:
			##data=data.copy()
			### ind=np.where(np.isfinite(data)==True)[0]
			### if vmin is None:
			###     vmin=np.min(data[ind])
			###     vmax=np.max(data[ind])
			###     vmin=vmin-(vmax-vmin)/255.0
			### else:
			###     ind1=ind[np.where(data[ind]<=vmin)[0]]
			###     ind2=np.where(np.isfinite(data)==False)[0]
			###     data[ind1]=vmin+(vmax-vmin)/255.0
			###     data[ind2]=vmin                
	
		##d=data.reshape(self.bins).transpose()
		##cmap=cm.get_cmap(cmapname)
		###cmap.set_bad('w',1.)
		###plt.imshow(np.ma.masked_where(np.isnan(d),d), interpolation=interpolation,extent=[self.range[0][0], self.range[0][1], self.range[1][0], self.range[1][1]],cmap=cmap,norm=norm,origin='lower',aspect='auto',vmin=vmin,vmax=vmax,**kwargs)
		##if smooth: 
			##kernel = astropy.convolution.Gaussian2DKernel(stddev=2.0)
			##d=astropy.convolution.convolve(d,kernel)
	###           d[d<1]=np.nan
			##kwargs['interpolation']='bicubic'
	
		##xe = np.linspace(self.range[0,0],self.range[0,1],self.bins[0]+1)
		##ye = np.linspace(self.range[1,0],self.range[1,1],self.bins[1]+1)
		##xc = (xe[0:-1]+xe[1:])/2.0
		##yc = (ye[0:-1]+ye[1:])/2.0
		##x=np.meshgrid(xc,yc,indexing='ij')
	
		##if ax is None:
			##ax=plt
		##if clevels is not None:
			##ax.contour(x[0],x[1],d.transpose(),clevels,vmin=vmin,vmax=vmax)
	
		##return ax.imshow(d,extent=[self.range[0][0], self.range[0][1], self.range[1][0], self.range[1][1]],cmap=cm.get_cmap(cmapname),origin='lower',aspect='auto',**kwargs)
	
	##def apply(self,quant,func,*args,**kwargs):
		##"""
		##Apply a scalar function to data points in a bin of the histogram.
		##Args:
	
			##quant: an array of size N (no of data points)
	
			##func: a scalar function to apply to data points          
	
		##Example:
		##>>> d_estimate=h.apply(d,np.mean)
			
		##"""
		##result=np.zeros(self.data.size,dtype=np.float64)+np.nan
		##for i in np.arange(self.data.size):
			##if self.data[i] > 0:
				##if 'weights' in kwargs:
					##result[i]=func(quant[self.indices(i)],*args,weights=kwargs['weights'][self.indices(i)])
				##else:
					##result[i]=func(quant[self.indices(i)],*args,**kwargs)
		##return result
				
	##def equisample(self,nsize):
		##indt=[]
		##for i in np.arange(self.data.size):
			##ind=self.indices(i)
			##np.random.shuffle(ind)
			##if ind.size>nsize:
				##indt.append(ind[0:nsize])
			##else:
				##indt.append(ind)
		##return np.concatenate(indt)
	
	##def resample(self,X,fsample=1.0,keys=None,verbose=False):
		##if type(X) == hist_nd:
			##hist1=X
		##else:
			##hist1=hist_nd(X,bins=self.bins,range=self.range,keys=keys,ri=True)
		##hist2=self
		##count1=np.sum(hist1.data)
		##count2=np.sum(hist2.data)
		##indices=np.zeros(np.int64(count2*fsample),dtype=np.int64)
		##nprev=0
		##misscount=0
		##for i in np.arange(hist1.data.size):
			##if hist2.data[i] > 0:
				##if hist1.data[i] < hist2.data[i]*fsample: 
					##misscount=misscount+hist2.data[i]*fsample-hist1.data[i]
				##if hist1.data[i] > 0:
					##ind1=hist1.indices(i)
					##np.random.shuffle(ind1)
					##nsize=np.int64(fsample*hist2.data[i])
					##ind=np.arange(nsize)
					##indices[nprev:nsize+nprev]=ind1[ind%ind1.size]
					##nprev=nsize+nprev
	
		##indices=indices[0:nprev]
		##if verbose:
			##print 'fsmaple=',fsample,' primary size=',np.sum(hist2.data),' secondary size=', np.sum(hist1.data)
			##print 'Points lacking in primary=',misscount
			##print 'Total points missed      =',int(np.sum(hist2.data)*fsample)-indices.size
		##return indices 
	
	
	##def compute(self,X,bins=10,range=None,normed=False,weights=None,keys=None,ri=True):
		##if type(keys) == str:
			##keys=[keys]
	
		##if type(X) == list:
			##dims=len(X)
			##V=X            
		##elif type(X) == tuple:
			##dims=len(X)
			##V=X            
		##elif type(X) == dict:
			##dims=len(keys)
			##V=[X[key] for key in keys]
		##elif type(X) == np.ndarray:
			##if X.dtype.char == 'V':
				##dims=len(keys)
				##V=[X[key] for key in keys]
			##else:
				##if X.ndim == 2:
					##dims=X.shape[1]
					##V=[X[:,i] for i in np.arange(dims)]
				##elif X.ndim == 1:
					##dims=1
					##V=[X]
				##else:
					##raise RuntimeError('Input data must be of ndim=1 or 2')
		
		##for i in np.arange(dims):
			##if type(V[i]) != np.ndarray:
				##V[i]=np.array(V[i])
	
		##nsize=V[0].size
		##for i in np.arange(dims):
			##if V[i].size != nsize:
				##print 'dim=',i,' size=',V[i].size,' required size=',nsize
				##raise RuntimeError('all dimensions not of same size')
		###print 'dims=',dims
	
		### V={}
		### if keys==None:
		###     if X.ndim == 2:
		###         dims=X.shape[1]
		###         for i in np.arange(dims):
		###             V[i]=X[:,i]
		###     elif X.ndim == 1:
		###         dims=1
		###         V[0]=X
		###     else:
		###         raise RuntimeError('Input data must be of ndim=1 or 2')
		### else:
		###     if type(keys) == str:
		###         keys=[keys]
		###     dims=len(keys)
		###     for i,key in enumerate(keys):
		###         V[i]=X[key]
	
		###determine range and bin size
		##flag=np.zeros(V[0].size,dtype=np.int64)
		##nbins=np.zeros(dims,dtype=np.int64)
		##nbins[:]=bins
		##total_bins=np.prod(nbins)
		##if range is None:
			##range=np.zeros((dims,2),dtype=np.float64)
			##for i in np.arange(dims):
				##range[i,0]=np.min(V[i])
				##range[i,1]=np.max(V[i])                
		##else:
			##range=np.float64(range)
			##if len(range.shape) == 1:
				##range=range.reshape((1,2))                    
			### clip the data to fit in range
			##for i in np.arange(dims):
				### To prevent last bin from taking extra points, on 22 April 2017 changed to from 
				### ind=np.where((V[i]<range[i,0])|(V[i]>range[i,1]))[0]
				##ind=np.where((V[i]<range[i,0])|(V[i]>=range[i,1]))[0]
				##flag[ind]=1                
		##bs=np.float64((range[:,1]-range[:,0]))/nbins
		##ind=np.where(flag==0)[0]
	
		###index of each data point
		##d=np.clip(np.int64((V[0]-range[0,0])/bs[0]),0,nbins[0]-1)
		##for i in np.arange(1,dims):                     
			### The scaled indices, 
			##d=nbins[i]*d+np.clip(np.int64((V[i]-range[i,0])/bs[i]),0,nbins[i]-1)
		##d=d[ind]
	
	###        h,e=np.histogram(d,range=[0,total_bins],bins=total_bins,normed=normed,weights=weights)
	
		### np.bincount is faster
		##h=np.bincount(d,weights=None,minlength=total_bins)
		##e=np.arange(total_bins+1)
	
		### reverse indices
		##if ri == True:
			##inds=np.argsort(d)
			##rev_ind=np.zeros(e.size+inds.size,dtype='int64') 
			##rev_ind[0]=e.size        
			##rev_ind[1:e.size]=np.cumsum(h)+e.size        
			##rev_ind[e.size:]=ind[inds]
		##else:
			##rev_ind=None
	
		##if weights is not None:
			##h=np.bincount(d,weights=weights[ind],minlength=total_bins)
	
		##if normed==True:
			##h=np.float64(h)/(np.sum(h)*np.prod(bs))
	
		##return h,nbins,range,rev_ind
	
	
	##def interp_weights(xyz, uvw):
	##d=xyz.shape[1]
	##tri = qhull.Delaunay(xyz)
	##simplex = tri.find_simplex(uvw)
	##vertices = np.take(tri.simplices, simplex, axis=0)
	##temp = np.take(tri.transform, simplex, axis=0)
	##delta = uvw - temp[:, d]
	##bary = np.einsum('njk,nk->nj', temp[:, :d, :], delta)
	##return vertices, np.hstack((bary, 1 - bary.sum(axis=1, keepdims=True)))
	
	##def interpolate(values, vtx, wts, fill_value=np.nan):
	##ret = np.einsum('nj,nj->n', np.take(values, vtx), wts)
	##if fill_value is not None:
		##ret[np.any(wts < 0, axis=1)] = fill_value
	##return ret
	
	##def npstruct_add_field(data1,fieldname,value):
	##if data1.dtype.char == 'V':
		##dt=[(name,data1.dtype[name]) for name in data1.dtype.names]
		##if fieldname not in data1.dtype.names:
			##dt.append((fieldname,value.dtype))
			##dt=np.dtype(dt)
			##data2=np.zeros(data1.size,dtype=dt)
			##for key in data1.dtype.names:
				##data2[key]=data1[key]
			##data2[fieldname]=value
			##return data2
		##else:
			##raise RuntimeError('Field already exists:'+fieldname)
			##return None
	##else:
		##raise RuntimeError('Data is not a numpy structure')
	
	
	##def delta(x):
	##return np.r_[x[1]-x[0],(x[2:]-x[0:-2])*0.5,x[-1]-x[-2]]
	
	
	##def ks(x1,x2,shift=None,tol=None):
	##ind=np.where(np.isfinite(x1))[0]
	##x1=x1[ind]
	##x1=np.sort(x1)
	##y1=np.arange(x1.size)*1.0/x1.size
	##if shift is not None:
		##xmin=np.min(x2)*np.min(shift)*0.99
		##xmax=np.max(x2)*np.max(shift)*1.01
	##else:
		##xmin=np.min(x2)*0.99
		##xmax=np.max(x2)*1.01
	##if xmin < np.min(x1):
		##x1=np.r_[xmin,x1]
		##y1=np.r_[0.0,y1]
	##if xmax> np.max(x1):
		##x1=np.r_[x1,xmax]
		##y1=np.r_[y1,1.0]        
	##func1=scipy.interpolate.interp1d(x1,y1)
	##ind=np.where(np.isfinite(x2))[0]
	##x2=x2[ind]
	##x2=np.sort(x2)    
	##y2=np.arange(x2.size)*1.0/x2.size
	##if shift is None:
		##return np.max(np.abs(func1(x2)-y2))*np.sqrt((x1.size*x2.size*1.0)/(x1.size+x2.size))
	##else:
		##gr=1.61803
		##r1=shift[0]
		##r4=shift[1]
		##r2=r4-(r4-r1)/gr
		##r3=r1+(r4-r1)/gr
		##s1=np.max(np.abs(func1(x2*r1)-y2))
		##s2=np.max(np.abs(func1(x2*r2)-y2))
		##s3=np.max(np.abs(func1(x2*r3)-y2))
		##s4=np.max(np.abs(func1(x2*r4)-y2))
		##while np.abs(r3-r2) > tol:
			##if s3 > s2:
				##r4=r3
				##s4=s3
				##r3=r2
				##s3=s2
				##r2=r4-(r4-r1)/gr
				##s2=np.max(np.abs(func1(x2*r2)-y2))
			##else:
				##r1=r2
				##s1=s2
				##r2=r3
				##s2=s3
				##r3=r1+(r4-r1)/gr
				##s3=np.max(np.abs(func1(x2*r3)-y2))                
		### res=[]
		### for temp in shift:        
		###     res.append(np.max(np.abs(func1(x2*temp)-y2)))
		##return r2,s2
	
	
	##def hist2d(x, y, bins=10, range=None, normed=False, weights=None, vmin=None, vmax=None, norm=None,cmapname='newjet'):
	##h1=hist_nd([x,y],range=range,normed=normed,weights=weights)
	##h1.imshow(cmapname=cmapname,vmin=vmin,vmax=vmax,norm=True)
	
	
	##class hist_ndp(hist_nd):    
	##def imshow1(self,density=False,vmin=None,vmax=None,**kwargs):
	###        putil.make_cmap('newjet')
		##data=self.data.copy()
		##if density:
			##bs=np.float64((self.range[:,1]-self.range[:,0]))/self.bins
			##data=np.float64(data)/(np.sum(data)*np.prod(bs))
		##else:
			##ind=np.where(np.isfinite(data)==True)[0]
			##if vmin is None:
				##vmin=np.min(data[ind])
				##vmax=np.max(data[ind])
				##vmin=vmin-(vmax-vmin)/255.0
			##else:
				##ind1=ind[np.where(data[ind]<=vmin)[0]]
				##ind2=np.where(np.isfinite(data)==False)[0]
				##data[ind1]=vmin+(vmax-vmin)/255.0
				##data[ind2]=vmin                
		##plt.imshow(data.reshape(self.bins).transpose(),extent=[self.range[0][0], self.range[0][1], self.range[1][0], self.range[1][1]],vmin=vmin,vmax=vmax,origin='lower',aspect='auto',**kwargs)
	### def recbin(x, y, C=None, gridsize=100, bins=None, xscale='linear', yscale='linear', extent=None, cmap=None, norm=None, vmin=None, vmax=None, alpha=None, linewidths=None, edgecolors='none', reduce_C_function=<function mean>, mincnt=None, marginals=False, hold=None, **kwargs):
	### h1=hist_nd([x,y],range=range,normed=normed,weights=weights)
	
	##class StructGrid():
	##"""
	##Takes an ndarray of npstruct and makes a grid out of it 
	##for use in interpolation.
	##"""
	##@staticmethod
	##def groupby(data,keys,regular=True):
		##"""
		##Will group data 
		##Parameters
		##----------
		##data :  (ndarray of npstruct)
		##keys : a list of strings specifying struct field names 
			   ##to group by
		##Returns
		##-------
		##res   : (list) data split in groups
		##groups: (list) of length len(keys)
				 ##groups[j][i] value of key[j] in i-th group         
		##"""
		##temp=0
		##for key in keys:
			##temp=temp+np.abs(np.diff(data[key]))
		##res=np.split(data,np.where(temp!=0)[0]+1)
	
		##groups1=[]
		##for key in keys:
			##groups1.append([])
		##for temp in res:
			##for i,key in enumerate(keys):
				##groups1[i].append(temp[key][0])
		##groups1=[np.array(group) for group in groups1]
	
		##if regular==False:
			##return res,groups1
		##else:
			##groups2=[np.unique(data[key]) for key in keys]
			##myshape=[group.size for group in groups2]
			##points=(np.meshgrid(*groups2,indexing='ij'))
			##groups2=[np.ravel(point) for point in points]
	
	
			##print 'groupby keys:',keys
			##print 'group size:',groups1[0].size
			##print 'grid size :',groups2[0].size,', ',myshape 
	
			##j=0
			##for i in range(groups2[0].size):
				##temp=0.0
				##for k in range(len(groups1)):
					##temp=temp+np.abs(groups1[k][j]-groups2[k][i])
				##if temp != 0:
					##res.insert(i,None)
				##else:
					##j=j+1
	
			##return [res,groups2]
	
	##@staticmethod
	##def fill_missing(data,keys,miss_order=None):        
		##"""
		
		##Parameters
		##----------
		##data :  (ndarray of npstruct)
		##keys : a list of strings specifying struct field names 
			   ##to group by
		##Returns
		##-------
		##res   : (list) data split in groups
		##groups: (list) of length len(keys)
				 ##groups[j][i] value of key[j] in i-th group         
		##"""
		##data.sort(order=keys)    
		##vnames=[temp for temp in data.dtype.names if temp not in keys]
		##myshape=np.array([np.unique(data[key]).size for key in keys])
		##ind=np.where(np.isfinite(data[vnames[0]])==False)[0]
		##if miss_order is None:
			##miss_order=range(len(keys))
		##else:
			##if type(miss_order[0]) == str:
				##miss_order=[keys.index(temp) for temp in miss_order]
		##print 'Filling missing grid points------>' 
		##print 'grid shape:',myshape
		##print 'missing ',ind.size,' out of ', data.size 
		##idone=0
		##idone1=0
		##for i in ind:
			##index=np.unravel_index(i,myshape)
			##ibk=None
			##for j in miss_order:
				##k=index[j]
				##i1=i
				##i2=i
				##while (k+1)<myshape[j]:
					##k=k+1
					##index1=list(index)
					##index1[j]=k
					##temp=np.ravel_multi_index(index1,myshape)
					##if np.isfinite(data[vnames[0]][temp]) == True:
						##i1=temp
						##break
				##k=index[j]
				##while (k-1)>0:
					##k=k-1
					##index1=list(index)
					##index1[j]=k
					##temp=np.ravel_multi_index(index1,myshape)
					##if np.isfinite(data[vnames[0]][temp]) == True:
						##i2=temp
						##break
				##if (i1!=i)and(i2!=i):
					##idone=idone+1
					##u=(data[keys[j]][i]-data[keys[j]][i1])/(data[keys[j]][i2]-data[keys[j]][i1])
					##for vname in vnames:
						##data[vname][i]=data[vname][i1]*(1-u)+data[vname][i2]*u
					##break
				##if ibk is None:
					##if i1!=i:
						##ibk=i1
					##elif i2 != i:
						##ibk=i2
			##if np.isfinite(data[vnames[0]][i]) == False:
				##if ibk is not None:
					##for vname in vnames:
						##data[vname][i]=data[vname][ibk]
					##idone=idone+1
					##idone1=idone1+1
	
		##print idone  ,' out of ', ind.size,' filled (', idone1, ' using nearest)'
	
	
	##@staticmethod
	##def griddata(data1,keys,qname,quant,miss_order=None):
		##"""
		##Construct a regular grid out of given points. 
		##One dimension can be regridded if it is not already regular.
		##If few grid points are missing it will also fill them. 
		##Parameters
		##----------
		##data1 : (ndarray of npstruct)
		##keys  : a list of strings specifying struct field names 
				##which will be grid dimensions
		##qname : field name which is irregular and needs to be regridded
				##qname should be present in list keys            
		##quant : the array specifying new values of qname
		##Returns
		##-------
		##data_new: (ndarray of npstruct) conforming to a regular grid,  
				  ##sorted and ordered by keys.
		##"""
		##data=np.resize(data1,data1.size)
		##if qname not in keys:
			##raise RuntimeError('qname should be present in keys')
		##split_keys=[temp for temp in keys if temp != qname]
		##vnames=[temp for temp in data.dtype.names if temp not in keys]
		##data.sort(order=split_keys+[qname])
		##data_split,groups=StructGrid.groupby(data,keys=split_keys)
	
		##mylist=[]
		##for j,d in enumerate(data_split):
			##data2=np.zeros(quant.size,dtype=data.dtype)
			##data2[qname]=quant
			##for i,key in enumerate(split_keys):
				##data2[key]=groups[i][j]
			##if d is None:
				##for vname in vnames:
					##data2[vname]=np.nan
			##else:
				##ind1=np.argsort(d[qname])                    
				##xs=np.r_[quant[0],d[qname][ind1],quant[-1]]
				##for i,vname in enumerate(vnames):
					##ys=np.r_[d[vname][ind1[0]],d[vname][ind1],d[vname][ind1[-1]]]
					##data2[vname]=scipy.interpolate.interp1d(xs,ys,assume_sorted=True)(data2[qname])
			##mylist.append(data2)
				
	
		##data_new=np.concatenate(mylist)
		##StructGrid.fill_missing(data_new,keys,miss_order)
		##data_new.sort(order=keys)
		##return data_new
	##class InterpGrid():
	##"""
	##Takes an ndarray of npstruct and makes a grid out of it 
	##for use in interpolation.
	##"""
	##def __init__(self,data1,keys):
		##"""
		##Parameters
		##----------
		##data1 :  (ndarray of npstruct)
		##keys : a list of strings specifying struct field names 
			   ##to group by
		##"""
		##data=np.resize(data1,data1.size)
		##data.sort(order=keys)
		##self.keys=keys
		##self.vnames=[temp for temp in data.dtype.names if temp not in self.keys]
	
		##self.points=[np.unique(data[key]) for key in self.keys]
		##self.delta=[self._get_delta(x) for x in self.points]
	
		##self.values={}        
		##for vname in self.vnames:
			##self.values[vname]=data[vname].reshape([point.size for point in self.points])
	
		##self.points1=tuple([data[key] for key in self.keys])
		##self.values1={}        
		##temp=np.meshgrid(*self.delta,indexing='ij')
		##self.vol=1.0
		##for x in temp:
			##self.vol=self.vol*x
		##for vname in self.vnames:
			##self.values1[vname]=data[vname]
	
	##def _get_delta(self,x):
		##return np.r_[x[1]-x[0],(x[0:-2]-x[2:])*0.5,x[-1]-x[-2]]
	
	##def _homogenize_arrays(self,xi):
		##xj=[np.asarray(t) for t in xi]
		##temp=xj[0]
		##for t in xj:
			##temp=temp+t
		##xj=[np.zeros_like(temp)+t for t in xj]
		##return xj
	
	##def get_values(self,vname,xi,fill_value=None,method='linear'):
		##"""
		##Parameters
		##----------
		##vname: field name whose interpolated value is desired
		##xi   : the coordinates for which interpolation is required
			   ##should be in same order as keys
		##Returns
		##-------
		##value: the interpolated value
		##"""
	
		##fill_value1=np.nan
		##if type(xi) == list:
			##xi=np.array(self._homogenize_arrays(xi)).transpose()
		##t1=scipy.interpolate.interpn(self.points,self.values[vname],xi,bounds_error=False,fill_value=fill_value1,method=method)
		##if fill_value == 'nearest':
			##ind=np.where(np.isfinite(t1)==False)[0]
			##if ind.size>0:
				##print 'outside interp range',ind.size,' out of ',t1.size 
                ##t1[ind]=scipy.interpolate.griddata(self.points1,self.values1[vname],xi[ind],method='nearest') 
        ##return t1
	
	
	##return 














