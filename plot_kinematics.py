## ==================================================================
## 		Plot the outputted dynamics maps
## ==================================================================
## warrenj 20150330 Process to plot the results of pPXF and GANDALF 
## routines.
## warrenj 20150727 Changing to a python script
## warrenj 20150917 Altered to plot and save all 8 plots.
## warrenj 20151216 Added section to plot residuals.
## warrenj 20160111 Add section to plot histgrams of the fields.
## warrenj 20160405 Added keyword CO to overlay CO maps from ALMA if avaible.
## This supersedes plot_results_CO.py
## warrenj 20160726 Now plots in a more object orientated way and creates a
## grid of plots too. This supersedes plot_results2.py

## *************************** KEYWORDS ************************* ##
# galaxy 		Name of the galaxy being plotted: used to find 
#				correct files and to print onto the plot.
# discard	0	Interger giving the number of rows and columns 
#				to be removed from the plot to remove edge 
#				effects.
# wav_range 	null	Imposed wavelength range on top of the automated 
#				limits.	
# vLimit 	2 	Integer giving the number of lowest and highest 
#				results in the plot to be discarded. Defualt 
#				ignores 2 highest and 2 lowest bins.
# norm		"lwv"	Normalisation methods for velocity fields:
#				lwv: luminosity weighted mean of the whole 
#				field is set to 0.
#				lum: velocity of the brightest spaxel is set 
#				to 0.
#				sig: Noralised to the mean velocity of 5 bins with the
#				highest LOSVD.
# plots 	False   Boolean to show plots as routine runs.
# nointerp 	False 	Boolean to use interpolation between bins in 
#				plots or not.
# residual 	False	Method to measure the residuals:
#			mean: use the mean of the residuals in each 
#				bin.
#			median: use the median of the residuals in 
#				each bin.
#			max: use the maximum of the residuals in 
#				each bin.
#			False: do not calculate and produce plot of 
#				residuals.
# CO	   False	Boolean to show ALMA CO plots overlaied (if they exist)
# D 		None Option to pass in the Data object instead of loading it.
## ************************************************************** ##

import numpy as np # for array handling
import glob # for searching for files
from astropy.io import fits # reads fits files (is from astropy)
from astropy import units as u
from astropy.coordinates import SkyCoord
from matplotlib.collections import LineCollection
import matplotlib.pyplot as plt # used for plotting
from plot_velfield_nointerp import plot_velfield_nointerp # for plotting with no interpolations. 
from plot_histogram import plot_histogram
import os
from sauron_colormap2 import sauron2 as sauron
import cPickle as pickle
from checkcomp import checkcomp
cc = checkcomp()

# # Give axes a saveTo property
# plt.axes.saveTo = property(lambda self:str())
# # Give axes an x and y on figure grid property
# plt.axes.figx = property(lambda self:int())
# plt.axes.figy = property(lambda self:int())
# # give axes a property to hold a colorbar axes
# plt.axes.cax = property(lambda self:plt.axes())
# # give axes a property to hold 2 additional axes for showing other axis
# plt.axes.ax2 = property(lambda self:plt.axes())
# plt.axes.ax3 = property(lambda self:plt.axes())

vin_dir = '%s/Data/vimos/analysis' % (cc.base_dir)
vin_dir_cube = '%s/Data/vimos/cubes' % (cc.base_dir)
ain_dir = '%s/Data/alma' % (cc.base_dir)
out_dir = '%s/Data/vimos/analysis' % (cc.base_dir)

#-----------------------------------------------------------------------------
def set_lims(v, positive=False, symmetric=False):
	if all(~np.isfinite(v)):
		return 0, 0

	v = v[np.isfinite(v)]

	for i in range(2):
		av = np.median(v)
		std = np.std(v)

		include = (v >= av - 3*std) * (v <= av + 3*std)
		v = v[include]

	vmin, vmax = min(v), max(v)

	if symmetric:
		vmax = np.mean([vmax, abs(vmin)])
		vmin = -vmax

	if positive:
		vmin = max(vmin, 0)

	return vmin, vmax
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
def add_CO(ax, galaxy, header, close=False):
	CO_image_dir="%s/%s-mom0.fits" % (ain_dir, galaxy)
	# Arcsec coords
	if os.path.exists(CO_image_dir) and not ax.RaDec:
		## *** This option doesn't seem to be working correctly and I don't understand 
		##	why... *******************************************************************
		alma = fits.open(CO_image_dir)[0]

		CO_x = np.arange(alma.header['NAXIS1'])*alma.header['CDELT1']*60*60
		CO_y = np.arange(alma.header['NAXIS2'])*alma.header['CDELT2']*60*60

		CO_x -= max(CO_x)/2
		CO_y -= max(CO_y)/2

		# Coordinates of VIMOS pointing
		vc = SkyCoord(header['HIERARCH CCD1 ESO INS IFU RA'], 
			header['HIERARCH CCD1 ESO INS IFU DEC'], 
			unit=(u.deg, u.deg))

		# Coordinates of ALMA pointing
		ac = SkyCoord(alma.header['CRVAL1'], alma.header['CRVAL2'],
			unit=(u.deg, u.deg))

		# Offset between the two pointings
		CO_x -= ((vc.ra.degree - header['CRPIX1']*header['CDELT1']/(60*60)) -
			(ac.ra.degree +
			alma.header['CRPIX1']*alma.header['CDELT1']/(60*60)))*60*60
				
		CO_y += ((vc.dec.degree - header['CRPIX2']*header['CDELT2']/(60*60)) -
			(ac.dec.degree +
			alma.header['CRPIX2']*alma.header['CDELT2']/(60*60)))*60*60

	#RA and dec coords
	elif os.path.exists(CO_image_dir):
		alma = fits.open(CO_image_dir)[0]

		ac = SkyCoord(alma.header['CRVAL1'], alma.header['CRVAL2'],
			unit=(u.deg, u.deg))

		CO_x = (np.arange(alma.header['NAXIS1']) - alma.header['CRPIX1']) *\
			alma.header['CDELT1'] + ac.ra.degree
		CO_y = (np.arange(alma.header['NAXIS2'])-alma.header['CRPIX2']) *\
			alma.header['CDELT2'] + ac.dec.degree

	# Plot and save
	if os.path.exists(CO_image_dir):
		#remove random extra dimenisons.
		CO_image = np.sum(np.sum(alma.data,axis=0), axis=0)
		cs = ax.contour(CO_x,CO_y,CO_image, colors='r')

		saveTo = os.path.dirname(ax.saveTo)+"/withCO/" + \
			os.path.basename(ax.saveTo)
		if not os.path.exists(os.path.dirname(saveTo)):
			os.makedirs(os.path.dirname(saveTo))
		plt.savefig(saveTo, bbox_inches="tight")

		if close:
			plt.close()
		else:
			# Make lines thinner for pdf by finding the line objects
			for o in ax.get_children():
				if type(o) is LineCollection:
					o.set_linewidth(0.3)
	plt.close('all')
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
def plot_results(galaxy, discard=0, wav_range="", vLimit=2, norm="lwv", 
	plots=False, residual=False, CO=False, show_bin_num=False,
	D=None, **kwargs):	

	SN_target = 30
	data_file =  "%s/galaxies.txt" % (vin_dir)
	# different data types need to be read separetly
	z_gals, SN_target_gals= np.loadtxt(data_file, unpack=True, skiprows=1, 
		usecols=(1,6))
	galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]
	z = z_gals[i_gal]
	SN_target=SN_target_gals[i_gal]


	if wav_range:
		wav_range_dir = wav_range + "/"
	else:
		wav_range_dir = ""

	# Return CO to False if ALMA file does not exist.
	if CO:
		CO_image_dir="%s/%s-mom0.fits" % (ain_dir, galaxy)
		if not os.path.exists(CO_image_dir): CO = False

	dataCubeDirectory = "%s/%s.cube.combined.corr.fits" % (vin_dir_cube, galaxy)
	output = "%s/%s/results/%s" % (out_dir, galaxy, wav_range_dir)
	out_nointerp = "%s/Documents/star_kin_paper/plots" % (cc.home_dir)
	# vin_dir_gasMC = "%s/%s/gas_MC" % (vin_dir, galaxy)
	out_pickle = '%s/pickled' % (output)

	# Used for CO plotting
	cubeFile = fits.open(dataCubeDirectory)
	header = cubeFile[0].header
	cubeFile.close()
# ------------== Reading pickle file and create plot  ===----------

	# Load pickle file from pickler.py
	if D is None:
		pickleFile = open("%s/dataObj_%s.pkl" % (out_pickle, wav_range), 'rb')
		D = pickle.load(pickleFile)
		pickleFile.close()

	# for bin in D.bin:
	# 	for l in bin.e_line.keys():
	# 		bin.components[l].__threshold__ = 3
	
	if D.norm_method != norm:
		D.norm_method = norm
		D.find_restFrame()

	
# ------------=============== Plot image ================----------
	
	print "    Image"
	
	title = "Total Flux"
	CBLabel = r"Flux (erg s$^{-1}$ cm$^{-2}$)"

	saveTo = "%s/%s_stellar_img.png" % (out_nointerp, galaxy)
	print saveTo
	fmin, fmax = set_lims(D.flux, positive=True)

	ax = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, D.flux, vmin=fmin, 
		vmax=fmax, nodots=True, show_bin_num=show_bin_num, colorbar=True, 
		label=CBLabel, cmap='gist_yarg', galaxy=galaxy.upper(), #header=header, 
		#signal_noise=D.SNRatio, signal_noise_target=SN_target, 
		save=saveTo)
	if CO:
		add_CO(ax, galaxy, header)

	if plots:
		plt.show()
# ------------========= Plot intensity (& EW) ===========----------
	print "    gas map(s) and equivalent widths"
	
	for c in D.e_components:
		print "        " + c

		if 'OIII' in c:
			c_title = '[OIII]'
		elif 'Hbeta' in c:
			c_title = r'H$_\beta$'
		elif 'Hgamma' in c:
			c_title = r'H$_\gamma$'
		else:
			c_title = c

		f_title = "%s Flux" % (c_title)
		# from header
		fCBtitle = r"Flux (erg s$^{-1}$ cm$^{-2}$)"
		f_min, f_max = set_lims(D.e_line[c].flux, positive=True)
		
		saveTo = "%s/%s_%s_img.png" % (out_nointerp, galaxy, c)
		
		ax = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, D.e_line[c].flux, 
			vmin=f_min, vmax=f_max, colorbar=True, nodots=True, label=fCBtitle, 
			cmap = 'gist_yarg', save=saveTo, galaxy=galaxy.upper()) #, header=header)
			
		if plots: plt.show()
		if CO:
			add_CO(ax, galaxy, header)
		

		eq_title = "%s Equivalent Width" % (c_title)
		eqCBtitle = r"Equivalent Width ($\AA$)"

		eq_min, eq_max = set_lims(D.e_line[c].equiv_width, positive=True)

		saveTo = "%s/%s_%s_eqW.png" % (out_nointerp, galaxy, c)

		ax = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, 
			D.e_line[c].equiv_width, vmin=eq_min, vmax=eq_max, colorbar=True, 
			nodots=True, label=eqCBtitle, galaxy=galaxy.upper(), #header=header, 
			save=saveTo)
		if CO:
			add_CO(ax, galaxy, header)
# ------------============ Amplitude/Noise ==============----------
		# amp_title = '%s Amplitude to Noise ratio' % (c_title)
		# amp_min, amp_max = set_lims(D.e_line[c].amp_noise, positive=True)
		# saveTo = "%s/%s_amp_nosie_%s.png" % (out_nointerp, c, wav_range)

		# ax1 = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, 
		# 	D.e_line[c].amp_noise, vmin=amp_min, vmax=amp_max, colorbar=True, 
		# 	nodots=True, title=amp_title, save=saveTo, close=not CO, header=header)
		# if CO:
		# 	ax1.saveTo = saveTo
		# 	add_CO(ax1, galaxy, header, close=True)
# ------------=========== Setting titles etc ============----------
	print '    Kinematics'
	# for c in ['stellar']: # For debugging
	for c in D.independent_components:
		print '        %s' % (c)
		im_type = c
		f = c
		if im_type == "gas":
			im_type=""
			f = 'Hbeta'
		elif im_type == "SF":
			im_type=" (Star Forming)"
			f = '[OIII]5007d'
		elif im_type == "Shocks":
			im_type=" (Shocking)"
			f = 'Hbeta'
		elif 'Hbeta' in im_type:
			im_type=" ("+r'H$_\beta$'+")"
		elif 'Hgamma' in im_type:
			im_type=" ("+r'H$_\gamma$'+")"
		elif 'OIII' in im_type:
			im_type=" (OIII)"
		else:
			im_type=" (" + im_type + ")"


		for k in D.components[f].plot.keys():
			symmetric=False
			positive=False
				
			CBLabel = None
			if k == "vel":
				title = 'Velocity'
				CBLabel = "V (km s$^{-1}$)"
				symmetric=True

			if  k == "sigma":
				title = 'Velocity Dispersion'
				CBLabel = r'$\mathrm{\sigma}$ (km s$^{-1}$)'
				positive = True

			if k == "h3":
				title = 'h3'
				symmetric = True

			if k == "h4":
				title = 'h4'


			if c == "stellar":
				utitle = "Stellar Uncertainty " + title + " Map"
				title = "Stellar " + title + " Map"
			else:
				utitle = "Ionised" + im_type + " Gas Uncertainty " + title + " Map"
				title = "Ionised" + im_type + " Gas\n" + title + " Map"
# ------------============ Setting v range ==============----------
			vmin, vmax = set_lims(D.components[f].plot[k], positive=positive, 
				symmetric=symmetric)
			v_uncert_min, v_uncert_max = set_lims(D.components[f].plot[k].uncert, 
				positive=True)

# ------------==== Plot velfield - no interperlation ====----------
			# Field plot
			saveTo = ("%s/%s_%s_%s.png" % (out_nointerp, galaxy, c, k))

			ax = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar,
				D.yBar, D.components[f].plot[k], vmin=vmin, vmax=vmax, #flux_type='notmag',
				nodots=True, show_bin_num=show_bin_num, colorbar=True, 
				label=CBLabel,galaxy = galaxy.upper(),
				signal_noise=D.SNRatio, #header=header, 
				signal_noise_target=SN_target, save=saveTo)
			if plots:
				plt.show()
			if CO:
				add_CO(ax, galaxy, header)
			
			# Uncertainty plot
			# saveTo = "%s/%s_%s_uncert_field_%s.png" % (out_nointerp, c, k, wav_range)
			# ax1 = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar,
			# 	D.components[c].plot[k].uncert, vmin=v_uncert_min, vmax=v_uncert_max,
			# 	flux_type='notmag', nodots=True, show_bin_num=show_bin_num,
			# 	colorbar=True, label=CBLabel, galaxy = galaxy.upper(),
			# 	redshift = z, title=utitle, save=saveTo, close=not CO, header=header)
			# if CO:
			# 	ax1.saveTo = saveTo
			# 	add_CO(ax1, galaxy, header, close=True)
				
			# if plots:
			# 	plt.show()
			#if CO:
			#	D.unbinned_flux = D.unbinned_flux_sav
# ------------============= Plot residuals ==============----------
	# if residual:
	# 	print "    " + residual + " residuals"

	# 	average_residuals = np.zeros(D.number_of_bins)
	# 	for i in range(D.number_of_bins):
	# 		bestfit = np.loadtxt('%s/bestfit/%d.dat' % (vin_dir_gasMC, i))
	# 		spectrum = np.loadtxt('%s/input/%d.dat' % (vin_dir_gasMC, i))
	# 		residuals = np.abs(spectrum - bestfit)
	# 		# remove edge pixels
	# 		residuals = np.delete(residuals, [np.arange(5), 
	# 		len(residuals)+np.arange(-5,0)], axis=0)

	# 		if residual=="mean":
	# 			average_residuals[i] = np.mean(residuals)
	# 		elif residual=="median":
	# 			average_residuals[i] = np.median(residuals)
	# 		elif residual=="max":
	# 			average_residuals[i] = np.max(np.abs(residuals))
				
	# 	minres, maxres = set_lims(average_residuals, positive=True) #mean_centered=True,
		
	# 	CBLabel = "Residuals"
	# 	title = str.capitalize(residual) + \
	# 	" Residuals of Bestfit to Normalised Spectrum"
	# 	saveTo = "%s/%s_residual_%s.png" % (out_nointerp, residual, wav_range)

	# 	ax1 = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar,
	# 		average_residuals, vmin=minres, vmax=maxres, flux_type='notmag',
	# 		nodots=True, show_bin_num=show_bin_num, colorbar=True, 
	# 		label=CBLabel, #flux_unbinned=D.unbinned_flux, 
	# 		galaxy = galaxy.upper(), redshift = z, title=title, 
	# 		save=saveTo, close=not CO, header=header)
	# 	if plots:
	# 		plt.show()
	# 	if CO:
	# 		ax1.saveTo = saveTo
	# 		add_CO(ax1, galaxy, header, close=True)
# # ------------=============== Plot Chi2/DOF =============----------
	# print "    chi2"

	# chi2 = np.zeros(D.number_of_bins)
	# for i in range(D.number_of_bins):
	# 	chi2[i] = np.loadtxt("%s/chi2/%d.dat" % (vin_dir_gasMC, i))

	# minchi2, maxchi2 = set_lims(chi2, positive = True)
	
	# CBLabel = "Chi2/DOF"
	# title = "Chi2/DOF of the bestfit"
	# saveTo = "%s/chi2_%s.png" % (out_nointerp, wav_range)

	# ax1 = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, chi2, 
	# 	vmin=minchi2, vmax=maxchi2, flux_type='notmag',
	# 	nodots=True, show_bin_num=show_bin_num, colorbar=True, 
	# 	label=CBLabel, flux_unbinned=D.unbinned_flux, 
	# 	galaxy = galaxy.upper(), redshift = z, title=title, 
	# 	save=saveTo, close=not CO, header=header)#, cmap=cm.blue)
	# if plots:
	# 	plt.show()
	# if CO:
	# 	ax1.saveTo = saveTo
	# 	add_CO(ax1, galaxy, header, close=True)
# ------------============ Line ratio maps ==============----------
	if any('OIII' in o for o in D.list_components):
		print "    line ratios"

		t_num = (len(D.e_components)-1)*len(D.e_components)/2
		for n in range(t_num):
			i = 0
			m = t_num
			while m > n:
				i += 1
				m -= i

			cA = D.e_components[len(D.e_components)-i-1]
			cB = D.e_components[len(D.e_components)-i+n-m]

			line_ratio = np.log10(D.e_line[cB].flux/D.e_line[cA].flux)
			if 'OIII' in cA:
				cA_title = '[OIII]'
			elif 'Hbeta' in cA:
				cA_title = r'H$_\beta$'
			elif 'Hdelta' in cA:
				cA_title = r'H$_\delta$'
			elif 'Hgamma' in cA:
				cA_title = r'H$_\gamma$'
			else:
				cA_title = cA

			if 'OIII' in cB:
				cB_title = '[OIII]'
			elif 'Hbeta' in cB:
				cB_title = r'H$_\beta$'
			elif 'Hdelta' in cB:
				cB_title = r'H$_\delta$'
			elif 'Hgamma' in cB:
				cB_title = r'H$_\gamma$'
			else:
				cB_title = cB
				
			lr_title = "%s/%s Line Ratio" % (cB_title, cA_title)
			lrCBtitle = r"log$_{10}$ (%s/%s)" %(cB_title,cA_title)

			lr_min, lr_max = set_lims(line_ratio)

			saveTo = "%s/lineratio/%s_%s_%s_LR.png" % (out_nointerp, galaxy, cB, cA)
			
			ax = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar,
				line_ratio, vmin=lr_min, vmax=lr_max, colorbar=True,
				nodots=True, label=lrCBtitle, #header=header, 
				show_bin_num=show_bin_num, galaxy = galaxy.upper(),
				save=saveTo)

			if CO:
				add_CO(ax, galaxy, header)



##############################################################################

# Use of plot_results.py

if __name__ == '__main__':

	galaxies = ['ic1459', 'ic1531', 'ic4296', 'ngc0612', 'ngc1399', 'ngc3100', 
		'ngc3557', 'ngc7075', 'pks0718-34', 'eso443-g024']
	# galaxies = [galaxies[8]]

	wav_range="4200-"
	discard = 2 # rows of pixels to discard- must have been the same 
			#	for all routines 
	vLimit = 2 #
	for galaxy in galaxies:
		print galaxy

		plot_results(galaxy, discard=discard, vLimit=vLimit, CO=True, norm='lwv',
			wav_range=wav_range, plots=False, residual = "median")