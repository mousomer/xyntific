# xyntific
**xyntific** is a fast lifting scheme impelementation of popular wavelet transforms (non dyadic DWT), written in low level C. This library is designed for non-dyadic data without resizing the input to dyadic length (i.e. powers of 2). Thus, non dyadic DWT can handle 39x39 images to maximum level, **without** intermediate buffering to 64x64.

Written and maintained by Omer Moussaffi (mousomer@gmail.com). Published under GPL 3.0.

Data files are assumed to be headerless, so user should retain lists of sizes (cols/rows) for data files. 

All wavelet transforms can be applied on any data size. Classical wavelet application is:

| Data:      | a0     | a1     | a2    | a3     | a4       | a5     | ...       | a(2n-3)   |a(2n-2)    |a(2n-1)    |
|:----------:|:------:|:------:|:-----:|:------:|:--------:|:------:|:---------:|:---------:|:---------:|:---------:|
| Transform  | low0   | high1  | low2  | high3  | low3     | ...    | low(2n-4) | high(2n-3)| low(2n-2) | high(2n-1)| 
| Pack       | low0   | low2   | low4  | ...    | low(2n-2)| high1  | high3     | high5     | ...       | high(2n-1)|

This implementation allows odd size data, and enables the transform to start from highpass. Perfect reconstruction is guaranteed as long as the user remembers to use the same boolean condition (startHigh = True/False) on the inverse transform. E.g.  

Transform with startHigh = false (5 lowpass elements, 4 highpass elements):

| Data:      | a0     | a1     | a2    | a3     | a4    | a5     | a6    | a7    |a8     |
|:----------:|:------:|:------:|:-----:|:------:|:-----:|:------:|:-----:|:-----:|:-----:|
| Transform  | low0   | high1  | low2  | high3  | low4  | high5  | low6  | high7 | low8  |
| Pack       | low0   | low2   | low4  | low6   | low8  | high1  | high3 | high5 | high7 |
| Inverse    | a0     | a1     | a2    | a3     | a4    | a5     | a6    | a7    |a8     |


or, Transform with startHigh = true (4 lowpass elements, 5 highpass elements):

| Data:      | a0     | a1     | a2    | a3     | a4    | a5    | a6    | a7    |a8     |
|:----------:|:------:|:------:|:-----:|:------:|:-----:|:-----:|:-----:|:-----:|:-----:|
| Transform  | high0  | low1   | high2 | low3   | high4 | low5  | high6 | low7  | high8 |
| Pack       | low1   | low3   | low5  | low7   | high0 | high2 | high4 | high6 | high8 |
| Inverse    | a0     | a1     | a2    | a3     | a4    | a5    | a6    | a7    |a8     |



Thus, for example, data of length 19 will have 10 highpass and 9 lowpass elements if startFroHigh is applied.

Implementation ispired by https://www.organicdesign.co.nz/Wavelet.c  
Command line parameter parsing borrowed from James Fowler's QccPack (http://qccpack.sourceforge.net/).  


# Structure
When building, include **upper library** --libxyn.h--

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


