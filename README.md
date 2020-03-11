# S&P 500 Volatility

This repository contains an exploratory R script related to past volatility in the S&P 500 as a predictor of future.

This is the hypothesis to test: the economic reaction to "news events" is incorporated into the price of the S&P 500 as observable sine waves in price, ("volatility,") decaying exponentially, as investors "buy low" and "sell high" in direct and indirect reaction to the news. An economically significant "news event" has a time of origin at which it "turns on," (representable by a Heaviside theta,) from which it oscillates into the future price of the index as a damped sine wave with node at the "turn on" point. News is either "good," (positive amplitude,) or "bad," (negative amplitude,) but either type of news results in relative oscillations of price between both higher and lower over time, as the ramifications of the news are priced into the index. These oscillations with origin times each have "frequencies," implying rates of time at which investors react to the news, but we do not hypothesize any model predicting the period on which investors react to any arbitrary economic news. So far, this work is exploratory, without quantitative validation.

Spectral analysis has been carried out on the S&P 500 price **differentials** at close-of-week since January 1st, 1928. Heaviside theta "turn-on" points have been omitted from the analysis so far, for model simplicity. So that we may linearly regress future price impacts, we incorporate the **exponential decay** of wave components into the regression model by way of an **exponential growth** transform on the S&P 500 index price, such that waves in volatility **never** decay in the transformed data, allowing for linear regression.

Hypothetically, there is a dependence of wave decay rate on wave frequency, akin to physical waves in a refractory medium. If the average APR of the S&P 500 price is (hypothetically, to basic general knowledge) a consistent 7% to 8%, then the lowest decay rate waves sustainable must have 7% to 8% **decay** rate. Lower decay rates would "inflate" exponentially faster than the average growth of the market, giving them imaginary number frequency and therefore only contributing to exponential growth or decay of price, therefore not exhibiting finite wave period. However, a **driven** wave with "driving" rate greater than 7% to 8% APR are also sustainable, for overcoming the "inflationary" drive of exponential market volume growth.

To sustain a wave, at the limit point of exactly infinite period for a driven wave, the exponential driving rate must limit to **just** infinite, (i.e., we assume the market is locally geometrically "flat," in the sense of Riemann). At the limit of the lifetime of the past market, **immediately**, a single lifetime past period drives the overall market average rate of return, (ex.: approximately 7% to 8% APR at the current epoch). At the limit of **0** frequency for market volatility, the **damping** coefficient must limit to exactly negatively infinite, or immediate total damping. This suggests that the damping rate for waves follows a form like d = c_0 * log(p / p_0), basically in analogy with the inverse of the continuously compounded interest formula, which we can render p = p_0 * e^(d/c_0) * e^t, where the damping rate is "d", the period of the volatility wave is "p," "c_0" is a proportionality constant to derive by inspection, and "p_0" is a characteristic volatility period.

We hope this repository will spur research and help give the market peace-of-mind!
