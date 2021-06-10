#include <stdlib.h>
#include <stdio.h>
#include <math.h>





double calculate_response(double freq, double sigma)
{
    /* magnitude response of discrete gaussian kernel       */
    /* in frequency domain                                  */
    /* input frequency is normalized so that the Nyquist    */
    /* frequency is pi                                      */

    double t = sigma * sigma;
    double response = exp(t*(cos(M_PI*freq)-1));

    return response;
}

double calculate_power_spectrum(double freq, double sigma)
{
    /* power spectrum is response magnitude squared         */
    /* input frequency is normalized so that the Nyquist    */
    /* frequency is pi                                      */

    double response = calculate_response(freq, sigma);

    return response * response;
}

double calculate_sigma(double cutoff_freq, double power_spectrum_level)
{
    /* calculate standard deviation of gaussian filter with */
    /* specified cutoff frequency. The cutoff frequency is  */
    /* defined as the frequency at which the power spectrum */
    /* drops below the specified level                      */
    /* input frequency is normalized so that the Nyquist    */
    /* frequency is pi                                      */
    /* the power spectrum level should be between o and 1   */


    double t = (0.5 * log(power_spectrum_level)) / (cos(M_PI * cutoff_freq) - 1.0);
    double sigma = sqrt(t);

    return sigma;
}

double calculate_cutoff(double sigma, double power_spectrum_level)
{
    /* calculate the cutoff frequency for a gaussian filter */
    /* with given standard deviation. The cutoff frequency  */ 
    /* is defined as the frequency at which the power       */
    /* spectrum drops below the specified level             */
    /* cutoff frequency is normalized so that the           */
    /* Nyquist frequency is pi                              */

    double t = sigma * sigma; 
    double cos_cutoff = (1/t) * (0.5 * log(power_spectrum_level)) + 1.0;
    double cutoff = acos(cos_cutoff) / M_PI;

    return cutoff;
}