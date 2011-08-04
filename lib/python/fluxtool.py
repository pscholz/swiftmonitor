#!/usr/bin/env python
# fluxtool written by Anne M. Archibald for the McGill Pulsar group, April 2007

from numpy import sqrt, sum, mean, shape, conjugate, argmax, where, amin, asarray, dot, amax, zeros, ones, floor, concatenate, arange, exp, pi, angle, all
from numpy.fft import rfft, irfft


class estimator(object):
    def __init__(self):
        pass

class rms_estimator(estimator):
    """RMS estimator of pulsed flux

    Estimate the pulsed flux in a pulse profile by smoothing it to a
    specified number of harmonics, then compute the RMS amplitude. 
    
    The normalization is chosen so that if smoothing is not a factor,
    the value returned is sqrt(mean(histogram**2)).
    """
    def __init__(self,n):
        estimator.__init__(self)
        self.name = "RMS estimator with %d harmonics" % n
        if n<1:
            raise ValueError, "Must include at least one harmonic"
        self.n = n
        self.RMS = True

    def __call__(self, histogram, uncertainties):

        if len(uncertainties)!=len(histogram):
            raise ValueError, "Uncertainties must be an array the same size as histogram"

        ft = rfft(histogram)

        # the variance on the real or imaginary part of a Fourier component
        # the 0th, and if n is even, the last, bin have zero imaginary part
        # so their variances are twice this.
        ft_variance_normal = sum(uncertainties**2)/2. 
        
        if self.n>=len(ft)-1:
            raise ValueError, "%d harmonics requested but only %d exist; this combination is not yet implemented" % (self.n, len(ft)-1)
        
        #power_estimate = sum(abs(ft[1:self.n+1])**2)
        power_estimate = max(sum(abs(ft[1:self.n+1])**2) - 2*self.n*ft_variance_normal,0)
        variance_estimate = ft_variance_normal

        normalization = len(histogram)**2/2.
        power_estimate /= normalization
        variance_estimate /= normalization

        return sqrt(power_estimate), sqrt(variance_estimate)

def cross_correlate(histogram, template, n_harmonics=None, upsample=16):
    """Find the shift required to align histogram with template

    """
    n = max(len(histogram),len(template))

    h_ft = rfft(histogram)
    t_ft = rfft(template)

    if len(h_ft)<len(t_ft):
        h_ft = concatenate((h_ft,zeros(len(t_ft)-len(h_ft))))
    elif len(t_ft)<len(h_ft):
        t_ft = concatenate((t_ft,zeros(len(h_ft)-len(t_ft))))

    if n_harmonics is not None:
        h_ft[n_harmonics+1:]*=0.
        t_ft[n_harmonics+1:]*=0.
    h_ft[0] = 0
    t_ft[0] = 0

    cross_correlations = irfft(conjugate(h_ft)*t_ft,n*upsample)
    shift = argmax(cross_correlations)/float(len(cross_correlations))

    assert 0<=shift<1

    #FIXME: warn if double-peaked

    return shift

def sum_u(values, uncertainties):
    values = asarray(values,dtype=float).ravel()
    uncertainties = asarray(uncertainties,dtype=float).ravel()
    if shape(values)!=shape(uncertainties):
        raise ValueError, "All arrays must be the same size"

    return (sum(values),
            sqrt(dot(uncertainties,uncertainties)))

def average(values, uncertainties, weights=None):
    """Compute the weighted average of a list of values with uncertainties

    Returns the average and its uncertainty.

    If weights are not specified, they are chosen so as to minimize 
    the resulting uncertainty.

    """
    values = asarray(values,dtype=float).ravel()
    uncertainties = asarray(uncertainties,dtype=float).ravel()
    if weights is None:
        weights = 1./uncertainties**2
    else:
        weights = asarray(weights,dtype=float).ravel()
    if shape(values)!=shape(uncertainties) or shape(values)!=shape(weights):
        raise ValueError, "All arrays must be the same size"

    return (dot(values,weights)/sum(weights),
            sqrt(dot(uncertainties**2,weights**2))/sum(weights))

class minimum_estimator(estimator):
    def __init__(self, template=None, off_pulse_bins=None, 
                       off_pulse_auto_margin=0., 
                       correlation_harmonics=None):

        self.template = template
        self.off_pulse_auto_margin = off_pulse_auto_margin

        if template is not None:
            self.name = "Minimum estimator using cross-correlator"
        elif off_pulse_bins is not None:
            self.name = "Minimum estimator using known minimum"
        else:
            self.name = "Minimum estimator finding minimum"
        self.RMS = False

        if off_pulse_bins is None:
            if template is not None:
                self.off_pulse_bins = list(where(
                    template-amin(template) <= off_pulse_auto_margin*(amax(template)-amin(template))
                    )[0])
            else:
                self.off_pulse_bins = None
        else:
            self.off_pulse_bins = list(off_pulse_bins)
            if not self.off_pulse_bins:
                raise ValueError, "No off-pulse bins specified!"
        self.correlation_harmonics = correlation_harmonics 

    def __call__(self, histogram, uncertainties):
        if self.template is not None:
            shift = cross_correlate(histogram, self.template, n_harmonics=self.correlation_harmonics)
        else:
            shift = 0.


        if self.off_pulse_bins is None:
            off_pulse_bins = list(where(
                histogram-amin(histogram) <= self.off_pulse_auto_margin*(amax(histogram)-amin(histogram))
                )[0])
        else:
            off_pulse_bins = self.off_pulse_bins

        if self.template is not None:
            temp_len = len(self.template)
        else:
            temp_len = len(histogram)
        # Okay, this sucks
        hist2 = concatenate((histogram,histogram))
        weights = zeros(len(hist2))
        for b in off_pulse_bins:
            if not 0<=b<temp_len:
                raise ValueError, "Invalid off-pulse bin %d" % b
            l = b/float(temp_len) - shift
            r = (b+1)/float(temp_len) - shift

            if l<0:
                k = floor(l)
                l-=k
                r-=k

            assert 0<=l
            assert l<r
            assert r-l<1.001/temp_len
            assert r<2

            for i in range(2*len(histogram)):
                lbe = i/float(len(histogram))
                rbe = (i+1)/float(len(histogram))
                if lbe<l:
                    if rbe>=l:
                        f = (rbe-l)/(rbe-lbe)
                        weights[i]+=f
                    else:
                        weights[i]+=0
                elif lbe<r:
                    if rbe>=r:
                        f = (r-lbe)/(rbe-lbe)
                        weights[i]+=f
                    else:
                        weights[i]+=1
                else: #lbe>=r
                    weights[i]+=0

        assert sum(weights)!=0

        u2 = concatenate((uncertainties,uncertainties))
        min, u_min = average(hist2,u2,weights)
        avg, u_avg = average(histogram,uncertainties,ones(len(histogram)))
        
        return avg-min, sqrt(u_min**2+u_avg**2)

class smoothed_minimum_estimator(minimum_estimator):
    def __init__(self, template=None, off_pulse_bins=None, 
                       off_pulse_auto_margin=0., 
                       harmonics=5,
                       upsample=64):

        self.template = template
        self.off_pulse_auto_margin = off_pulse_auto_margin

        if template is not None:
            self.name = "Smoothed minimum estimator using cross-correlator"
        elif off_pulse_bins is not None:
            self.name = "Smoothed minimum estimator using known minimum"
        else:
            self.name = "Smoothed minimum estimator finding minimum"
        self.RMS = False

        self.harmonics = harmonics
        self.upsample = upsample
        self.name += " (%d harmonics)" % harmonics

        if off_pulse_bins is None:
            self.off_pulse_bins = None
        else:
            self.off_pulse_bins = list(off_pulse_bins)
            if not self.off_pulse_bins:
                raise ValueError, "No off-pulse bins specified!"

        if template is not None:
            self.off_pulse_ranges = self.compute_off_pulse_ranges(template,len(template))
        else:
            self.off_pulse_ranges = None

    def compute_off_pulse_ranges(self,histogram,n,peak=False):
        n = float(n)
        if self.off_pulse_bins is not None:
            if peak:
                off_pulse_bins = self.peak_bins
            else:
                off_pulse_bins = self.off_pulse_bins
        else:
            if peak:
                histogram = -histogram
            off_pulse_bins = list(where(
                    histogram-amin(histogram) <= 
                       self.off_pulse_auto_margin*(amax(histogram)-amin(histogram))
                    )[0])
        #if [a for a in off_pulse_bins if a>=n]:
        #    raise ValueError, "Specified bins '%s' include some larger than %d" % (off_pulse_bins, n)
        return [(b/n,(b+1)/n) for b in off_pulse_bins]
        
    def evaluate(self, ft, ft_variance_normal, lh, ranges):
        
        # this weird voodoo allows us to calculate not only the integral
        # of the smoothed function over the off-pulse interval but error
        # bars on the result you get.
        exponential_integral = zeros(self.harmonics,complex)
        width = 0.
        for (l,r) in ranges:
            k = arange(1,self.harmonics+1)
            exponential_integral += (exp(2.j*pi*k*r)-exp(2.j*pi*k*l))/(2.j*pi*k)
            width += r-l

        assert width>0

        # len(histogram) is there to deal with the normalization of numpy's rfft
        normalization = 2./lh
        min = normalization*dot(ft[1:self.harmonics+1],exponential_integral).real/width
        u_min = normalization*sqrt(sum(ft_variance_normal*abs(exponential_integral)**2))/width
        
        return min, u_min

    def find_shift(self, histogram):
        if self.template is not None:
            shift = cross_correlate(histogram, self.template, n_harmonics=self.harmonics)
            #print "shift:", shift
        else:
            shift = 0.
        return shift

    def __call__(self, histogram, uncertainties):
        shift = self.find_shift(histogram)

        ft = rfft(histogram)
        ft[self.harmonics+1:]*=0.
        ft_variance_normal = sum(uncertainties**2)/2. 

        if self.off_pulse_ranges is not None:
            off_pulse_ranges = [(l-shift, r-shift) for (l,r) in self.off_pulse_ranges]
        else: 
            off_pulse_ranges = self.compute_off_pulse_ranges(irfft(ft,self.harmonics*16),len(histogram))

        min, u_min = self.evaluate(ft,ft_variance_normal,len(histogram),off_pulse_ranges)
        # We left out the 0th Fourier coefficient when computing min, so we have essentially forced the average to be zero
        return -min, u_min
        #avg, u_avg = average(histogram,uncertainties,ones(len(histogram)))
        #return avg-min, sqrt(u_min**2+u_avg**2)


class peak_to_peak_estimator(smoothed_minimum_estimator):

    def __init__(self, template=None, off_pulse_bins=None, peak_bins=None,
                       off_pulse_auto_margin=0., 
                       harmonics=5,
                       upsample=64):

        self.template = template
        self.off_pulse_auto_margin = off_pulse_auto_margin

        if template is not None:
            self.name = "Smoothed peak-to-peak estimator using cross-correlator"
        elif off_pulse_bins is not None:
            self.name = "Smoothed peak-to-peak estimator using known minimum"
        else:
            self.name = "Smoothed peak-to-peak estimator finding minimum"
        self.RMS = False

        self.harmonics = harmonics
        self.upsample = upsample
        self.name += " (%d harmonics)" % harmonics

        if off_pulse_bins is None:
            if peak_bins is not None:
                raise ValueError, "Peak bins specified but off-pulse bins not specified!"
            self.off_pulse_bins = None
        else:
            self.off_pulse_bins = list(off_pulse_bins)
            if not self.off_pulse_bins:
                raise ValueError, "No off-pulse bins specified!"
        if peak_bins is None:
            if off_pulse_bins is not None:
                raise ValueError, "Off-pulse bins specified but peak bins not specified!"
            self.peak_bins = None
        else:
            self.peak_bins = list(peak_bins)
            if not self.peak_bins:
                raise ValueError, "No peak bins specified!"

        if template is not None:
            self.off_pulse_ranges = self.compute_off_pulse_ranges(template,len(template))
            self.peak_ranges = self.compute_off_pulse_ranges(template,len(template),peak=True)
        else:
            self.off_pulse_ranges = None
            self.peak_ranges = None


    def __call__(self,histogram,uncertainties):
        shift = self.find_shift(histogram)

        ft = rfft(histogram)
        ft[self.harmonics+1:]*=0.
        ft_variance_normal = sum(uncertainties**2)/2. 


        if self.off_pulse_ranges is None:
            off_pulse_ranges = self.compute_off_pulse_ranges(irfft(ft,self.harmonics*16),len(histogram))
        else:
            off_pulse_ranges = [(l-shift, r-shift) for (l,r) in self.off_pulse_ranges]
        
        if self.peak_ranges is None:
            peak_ranges = self.compute_off_pulse_ranges(irfft(ft,self.harmonics*16),len(histogram),peak=True)
        else:
            peak_ranges = [(l-shift, r-shift) for (l,r) in self.peak_ranges]

        min, u_min = self.evaluate(ft,ft_variance_normal,len(histogram),off_pulse_ranges)
        max, u_max = self.evaluate(ft,ft_variance_normal,len(histogram),peak_ranges)
        #raise ValueError((abs(ft[1]), angle(ft[1]), min, max, off_pulse_ranges, peak_ranges ))
        # We left out the 0th Fourier coefficient when computing min, so we have essentially forced the average to be zero
        return max-min, sqrt(u_max**2+u_min**2)
        #avg, u_avg = average(histogram,uncertainties,ones(len(histogram)))
        #return avg-min, sqrt(u_min**2+u_avg**2)


if __name__=='__main__':
    import sys
    from optparse import OptionParser
    import scipy.io
    parser = OptionParser("Usage: %prog [options] input_profile", version="%prog 1.0")
    parser.add_option("-n", "--harmonics", 
                      dest="harmonics", type="int", default=5,
                      help="Smooth to N harmonics (default %default)", 
                      metavar="N")
    parser.add_option("-t", "--template",
                      dest="template",
                      help="Load FILE for use as a template for cross-correlation",
                      metavar="FILE")
    parser.add_option("-f", "--off-pulse-bins",
                      dest="off_pulse_bins",
                      help="Treat LIST as a (comma-separated) list of bins in the template to treat as off-pulse time (may include ranges, as 12,15-17)",
                      metavar="LIST")
    parser.add_option("-a", "--automatic-off-pulse",
                      dest="auto_off_pulse", action="store_true",
                      help="Treat the lowest bin(s) of the template as off-pulse time")
    parser.add_option("-A", "--automatic-off-pulse-threshold",
                      dest="auto_off_pulse_threshold", type="float", default=0,
                      help="Treat as off-pulse bins any bin in the template that is within T*(max-min) of the minimum bin (default %default)", 
                      metavar="T")
    parser.add_option("-N", "--no-smooth-minimum", 
                      action="store_false", dest="smooth",
                      help="Turn off profile smoothing",
                      default=True)
    parser.add_option("-T", "--two-cycles", 
                      action="store_true", dest="twocycles",
                      help="Input file contains two cycles",
                      default=False)

    (options, args) = parser.parse_args()

    if len(args)!=1:
        parser.error("Exactly one positional argument required")
                      
    infile, = args

    if infile=="-":
        infile = sys.stdin

    try:
        profile = scipy.io.read_array(infile)
        n, m = shape(profile)
        if m!=3:
            raise ValueError
    except IOError:
        parser.error("Unable to read input file %s" % infile)
    except ValueError:
        parser.error("Input file must have three columns: bin number, count rate, and uncertainty")


    if not options.twocycles:
        histogram = profile[:,1]
        uncertainties = profile[:,2]
        if (len(histogram)%2==0 and 
            all(histogram[:len(histogram)//2]==
                   histogram[len(histogram)//2:])):
            sys.stderr.write("Warning: profile appears to contain two cycles\n")
    else:
        histogram = profile[:,1]
        uncertainties = profile[:,2]
        if len(histogram)%2==1:
            parser.error("Profile was supposed to contain two cycles but has odd length")
        if not all(histogram[:len(histogram)//2]==
                   histogram[len(histogram)//2:]):
            sys.stderr.write("Warning: profile does not appear to contain two cycles\n")
        uncertainties = uncertainties[:len(histogram)//2]    
        histogram = histogram[:len(histogram)//2]    


    total_flux = mean(histogram)

    rms_value, rms_uncertainty = rms_estimator(options.harmonics)(histogram, uncertainties)

    print "RMS pulsed flux:            \t%#0.7g\t+/-\t%#0.7g" % (rms_value, rms_uncertainty)
    print "RMS pulsed fraction:        \t%#0.7g\t+/-\t%#0.7g" % (rms_value/total_flux, rms_uncertainty/total_flux)

    compute_area = False
    if options.template is not None:
        compute_area = True
        try:
            profile = scipy.io.read_array(options.template)
            n, m = shape(profile)
            if m not in [2,3]:
                raise ValueError
        except IOError:
            parser.error("Unable to read template file %s" % infile)
        except ValueError:
            parser.error("Template file must have two or three columns: bin number, count rate, and optionally uncertainty")
        template = profile[:,1]

        t_value, t_uncertainty = rms_estimator(options.harmonics)(template, zeros(len(template)))
        scal = mean(template-amin(template))/t_value
        print "RMS to area scaling factor: \t\t%#0.7g" % scal
        print "Scaled RMS pulsed flux:     \t%#0.7g\t+/-\t%#0.7g" % (rms_value*scal, rms_uncertainty*scal)
        print "Scaled RMS pulsed fraction: \t%#0.7g\t+/-\t%#0.7g" % (rms_value*scal/total_flux, rms_uncertainty*scal/total_flux)

        
    else:
        template = None
    
    off_pulse_bins = None
    off_pulse_threshold = None

    if options.off_pulse_bins is not None: 
        compute_area = True
        def convert_bins(s):
            if "-" in s:
                a,b = s.split("-")
                return range(int(a),int(b)+1)
            else:
                return int(s),
        off_pulse_bins = concatenate(map(convert_bins,options.off_pulse_bins.split(",")))
        if options.auto_off_pulse:
            parser.error("Automatic off-pulse calculation is incompatible with an explicit list of off-pulse bins")
    elif options.auto_off_pulse:
        compute_area = True
        off_pulse_threshold = options.auto_off_pulse_threshold

    if compute_area and off_pulse_bins is None and off_pulse_threshold is None:
        parser.error("No method specified to choose off-pulse bins, please use either --off-pulse-bins or --automatic-off-pulse")
        compute_area = False

    if compute_area:    
        if options.smooth:
            E = smoothed_minimum_estimator(template,
                                           off_pulse_bins,
                                           off_pulse_threshold,
                                           options.harmonics)
        else:
            E = minimum_estimator(template,
                                  off_pulse_bins,
                                  off_pulse_threshold,
                                  None)
        (v, u) = E(histogram, uncertainties)
        print "Area pulsed flux:           \t%#0.7g\t+/-\t%#0.7g" % (v,u)
        print "Area pulsed fraction:       \t%#0.7g\t+/-\t%#0.7g" % (v/total_flux,u/total_flux)
