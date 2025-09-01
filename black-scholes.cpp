#include <iostream>
#include <cmath>

/**
 * @struct Contract
 * @brief Represents an options contract in the context of the Black-Scholes model.
 * 
 * @var premium
 * The theoretical price of the option calculated using the Black-Scholes formula.
 * 
 * @var dte
 * Days till expiry (DTE) of the option, calculated as the time to maturity (T) in years multiplied by 365.2425.
 * 
 * @var delta
 * The rate of change of the option's price with respect to changes in the price of the underlying asset.
 * 
 * @var gamma
 * The rate of change of delta with respect to changes in the price of the underlying asset.
 * 
 * @var theta
 * The rate of change of the option's price with respect to the passage of time (time decay).
 * 
 * @var vega
 * The rate of change of the option's price with respect to changes in the implied volatility of the underlying asset.
 * 
 * @var rho
 * The rate of change of the option's price with respect to changes in the risk-free interest rate.
 * 
 * @var implied_volatility
 * The volatility of the underlying asset implied by the market price of the option, estimated iteratively.
 * 
 * @var intrinsic_value
 * The value of the option if it were exercised immediately, calculated as the difference between the underlying asset price and the strike price (for calls) or vice versa (for puts), or zero if the option is out of the money.
 */
struct Contract
{
    double premium;
    int dte;
    double delta;
    double gamma;
    double theta;
    double vega;
    double rho;
    double implied_volatility;
    double intrinsic_value;
};


// Error function approximation for the cumulative standard normal density function
double erf(double x){
    const double A1 = 0.254829592;
    const double A2 = -0.284496736;
    const double A3 = 1.421413741;
    const double A4 = -1.453152027;
    const double A5 = 1.061405429;
    const double P = 0.3275911;

    // Save the sign of x
    int sign = (x >= 0) ? 1 : -1;
    x = fabs(x);

    // A&S formula 7.1.26
    double t = 1.0 / (1.0 + P * x);
    double y = 1.0 - (((((A5 * t + A4) * t) + A3) * t + A2) * t + A1) * t * exp(-x * x);

    return sign * y;
}

// Cumulative standard normal density function
double cumulativeStandardNormal(double x) {
    return 0.5 * (1.0 + erf(x / sqrt(2.0)));
}


// Standard normal PDF
double standardNormalPDF(double x) {
    return (1.0 / std::sqrt(2 * M_PI)) * std::exp(-0.5 * x * x);
}

// Black-Scholes pricing with embedded implied volatility calculation
Contract blackScholesOptionPricing(double S0, double K, double r, double sigma, double T, bool isCallOption) {
    Contract con;
    con.dte = static_cast<int>(T * 365.2425);

    double d1 = (log(S0 / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    double d2 = d1 - sigma * std::sqrt(T);

    if (isCallOption) {
        con.premium = S0 * cumulativeStandardNormal(d1) - K * std::exp(-r * T) * cumulativeStandardNormal(d2);
        con.delta = cumulativeStandardNormal(d1);
        con.gamma = standardNormalPDF(d1) / (S0 * sigma * std::sqrt(T));
        con.theta = -(S0 * standardNormalPDF(d1) * sigma) / (2 * std::sqrt(T)) - r * K * std::exp(-r * T) * cumulativeStandardNormal(d2);
        con.vega = S0 * standardNormalPDF(d1) * std::sqrt(T);
        con.rho = K * T * std::exp(-r * T) * cumulativeStandardNormal(d2);
        con.intrinsic_value = std::max(S0 - K, 0.0);
    } else {
        con.premium = K * std::exp(-r * T) * cumulativeStandardNormal(-d2) - S0 * cumulativeStandardNormal(-d1);
        con.delta = cumulativeStandardNormal(d1) - 1;
        con.gamma = standardNormalPDF(d1) / (S0 * sigma * std::sqrt(T));
        con.theta = -(S0 * standardNormalPDF(d1) * sigma) / (2 * std::sqrt(T)) + r * K * std::exp(-r * T) * cumulativeStandardNormal(-d2);
        con.vega = S0 * standardNormalPDF(d1) * std::sqrt(T);
        con.rho = -K * T * std::exp(-r * T) * cumulativeStandardNormal(-d2);
        con.intrinsic_value = std::max(K - S0, 0.0);
    }

    double market_price = con.premium;
    double tol = 1e-6;
    double sigma_guess = sigma;
    double price, vega, diff;

    // Calculate Implied Volatility
    do {
        double d1_nr = (log(S0 / K) + (r + 0.5 * sigma_guess * sigma_guess) * T) / (sigma_guess * std::sqrt(T));
        double d2_nr = d1_nr - sigma_guess * std::sqrt(T);

        if (isCallOption)
            price = S0 * cumulativeStandardNormal(d1_nr) - K * std::exp(-r * T) * cumulativeStandardNormal(d2_nr);
        else
            price = K * std::exp(-r * T) * cumulativeStandardNormal(-d2_nr) - S0 * cumulativeStandardNormal(-d1_nr);

        vega = S0 * standardNormalPDF(d1_nr) * std::sqrt(T);

        diff = market_price - price;
        sigma_guess += diff / vega;
    } while (fabs(diff) > tol);

    con.implied_volatility = sigma_guess;

    return con;
}


// Example
int main() {
    double S0 = 100.0;   // Initial stock price
    double K = 100.0;    // Strike price
    double r = 0.05;     // Risk-free rate
    double sigma = 0.2;  // Volatility
    double T = 1;      // Time to maturity (in years)
    

    // Calculate option prices
    auto callContract = blackScholesOptionPricing(S0, K, r, sigma, T, true);
    auto putContract = blackScholesOptionPricing(S0, K, r, sigma, T, false);

    // Output the results
    std::cout << "European Call Option Price: " << callContract.premium << "\ndte: " << callContract.dte << "\ndelta: " << callContract.delta << "\ngamma: " << callContract.gamma << "\ntheta: " << callContract.theta << "\nvega: " << callContract.vega  << "\nrho: " << callContract.rho << "\nimplied volatility: " << callContract.implied_volatility  << "\nintrinsic value: " << callContract.intrinsic_value << std::endl;
    std::cout << std::endl;
    std::cout << "European Put Option Price: " << putContract.premium << "\ndte: " << putContract.dte << "\ndelta: " << putContract.delta << "\ngamma: " << putContract.gamma << "\ntheta: " << putContract.theta << "\nvega: " << putContract.vega << "\nrho: " << putContract.rho << "\nimplied volatility: " << putContract.implied_volatility  << "\nintrinsic value: " << putContract.intrinsic_value << std::endl;
 
    return 0;
}