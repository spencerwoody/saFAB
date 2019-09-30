##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Efron
##' @param z 
##' @param nmids 
##' @param pct 
##' @param pct0 
##' @param df 
##' @param nulltype 
##' @return 
##' @author Spencer Woody
##'
##' @export
efron = function(z, nmids=150, pct=-0.01, pct0=0.25, df=10, nulltype='theoretical') {
	# estimate f(z) and f_0(z) using Efron (2004)'s method
	stopifnot(any(nulltype == 'theoretical', nulltype=='empirical'))
	
	result = tryCatch({

		N = length(z)
		med = median(z)
		myrange = med + (1 - pct) * (range(z) - med)
		lb = myrange[1]
		ub = myrange[2]
		                
		breaks = seq(lb, ub, length= nmids +1)
		h1 = hist(z, breaks = breaks, plot = FALSE)
		mids = (breaks[-1] + breaks[-(nmids+1)])/2
		zcounts = h1$counts
		glm1 = glm(zcounts ~ splines::ns(mids, df = df), family=poisson)
		zrate = glm1$fit
		D = (zcounts - zrate)/sqrt(zrate+1)
		D = D[-c(1,nmids)]
		if (sum(D^2) > qchisq(0.9, nmids-2-df)) {
			warning(paste0("f(z) misfit = ", round(D, 1), ".  Rerun with increased df."))
		}	
		
		zdens = {zrate/sum(zrate)}/diff(breaks)

		# Now do spline interpolation for the density at the observed points
		ispl2 = splines::interpSpline( zdens ~ mids )
		fz = predict(ispl2, z)$y
		
		# Pick out the middle of the data points
		ql = quantile(z, pct0)
		qu = quantile(z, 1-pct0)
		ind0 = intersect(which(z > ql), which(z<qu))
		if(nulltype=='empirical') {
		# empirical null by central moment matching
			z0 = z[ind0]
			l0 = log(fz[ind0])
			zmax = z[which.max(l0)]
			lm0 = lm(l0~I(z0-zmax) + I((z0-zmax)^2))
			b0 = coef(lm0)
			sig = as.numeric(sqrt(-1.0/{2.0*b0[3]}))
			mu = as.numeric(-b0[2]/(2.0*b0[3]) + zmax)
			# lm0 = lm(l0 ~ z0 + I(z0^2))
			# b0 = coef(lm0)
			# sig = as.numeric(sqrt(-1/{2*b0[3]}))
			# mu = as.numeric(b0[2] * sig^2)
		} else {
		# theoretical null
			sig = 1
			mu = 0
		}
		p0 = sum(fz[ind0])/sum(dnorm(z[ind0], mu, sig))
		localfdr = pmin(1, p0*dnorm(z, mu, sig)/fz)
		list(mids=mids, breaks=breaks, zcounts=zcounts, zdens=zdens,
			z=z, fz=fz, mu0=mu, sig0=sig, p0=p0, fdr=localfdr)
	}, error = function(err) {
		print(err)
		list(mids=NULL, breaks=NULL, zcounts=NULL, zdens=NULL,
			z=NULL, fz=NULL, mu0=NULL, sig0=NULL, p0=NULL, fdr=NULL)
	}, finally = {
		# Nothing to clean up
	})
	return(result)
}
