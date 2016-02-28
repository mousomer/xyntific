# xyntific
**xyntific** is a fast lifting scheme impelementation of popular wavelet transforms, written in low level C. This library is designed for non-dyadic data IO (e.g. can handle 39x39 images).

Data files are assumed to be headerless, so user should retain lists of sizes (cols/rows) for data files. 

All wavelet transforms can be applied on any data size. Classical wavelet application is:

Data a0 a1 a2 a3 a4 a5....

--> low0 high1 low2 high3....

This implementation can allow to start from highpass.

Thus, for example, data of length 9 will have 5 highpass and 4 lowpass elements if startFroHigh is applied.

Implementation uses some tools from James Fowler's QccPack (http://qccpack.sourceforge.net/).  
Ispired by https://www.organicdesign.co.nz/Wavelet.c

# Structure
**upper library**: --libxyn.h--

    general options for all executables:  
        -v verbose - display verbose information (mainly I/O)
        -type d|f|l|g|c - choose data type (integer|float|long|double exact|char). Currently, most functions support only d (integer), g (double) and f (float).  
        -length - integer length of data
        -fileName - input/output file name string  


    genRandFile.c [-v ] [-type dataType] [-min min] [-max max] length fileName  
        generates a ranadom number vector, between $min and $max, of length $length and puts it into a file name $fileName.  


    genLinearFile [-v ] [-type dataType] [-power power] [-min min] max len fileName  
        generates a polynomial file - y=a*x^n. $power is n (defualt 1), min and max are image min/max (def min=0)
        use these to test perfect reconstruction and statistical moments


    calcMoment [-v ] [-type dataType] [-central ] [-divImp ] [-length len] moment inpName  
        calculate the n-th statistical moment on the data.  
        -cenrtal for central moment (remove mean from data)  
        -divImp calls an implementation dividing data inplace (diffenent rounding errors)  
        -length inputs only $len of the data  
        moment is the statistical moment order (=n).  
        non-central output is double = sum(data^n)/n. cantral output is sum((data-mean)^n)/n.  
        mean energy is non-central second moment.  


    displayArray [-prec precision] [-digits digits] [-line lineLen] [-type  dataType(d|c|l|f|e)] [-length  len] fName
        displays the numeric contents of a binary file. 
        -prec is the number of dispayed fractional digits. Only relevant to floating point.
        -digits sets the spacde reserved for each number
        -line sets amount of numbers displayed per line.


    displayWavelet [-lift ] waveName
        shows wavelet coefficients.
        -lift shows lifting coefficients
        waveNames supported: Haar, CDF-5-3, CDF-9-7


    dwt1d/idwt1d [-v ] [-out outName] [-noscl ] [-type dataType] [-lift ] [-inPlace ] [-startHigh ] [-length len] waveName inpName
        -out output file name (default: inpName.fwd)
        -noscl - do not upscale integer data
        -lift call lifting implementation (usually better)
        -inPlace do not shuffle lowpass/highpass
        -startHigh - set first element of output to be highpass


    dwt2d/idwt2d [-v ] [-out outName] [-type dataType] [-lift ] [-inPlace ] [-startHigh ] [-nC nC] nR waveName inpName
        seperable 2d wavelet transform - preform 1d transform and rows and cols
        -nC sets number of columns (default is calculated from filesize/nR)
        nR number of data rows


