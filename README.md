# Repository for dengue trial project - MMED 2023

We aim to estimate the required sample size for an intervention trial in Gampaha, Sri Lanka (SLaDe trial).
The required sample size depends on the number of outcome events that are recorded.
It is most efficient to record outcomes that occur most frequently, so recording infections is more efficient than recording clinical cases.
But we do not know how many people get infected each year in Gampaha.
We only know the number of reported dengue cases in Gampaha (but not the % of infections that lead to clinical cases).

However, in Colombo we know both the annual probability of first infection and the number of reported dengue cases.
We use this to estimate the number of infections in Gampaha.
First, we calculate the FOI from the prob of first infection. We use this to estimate the number of people experiencing their first, second, third and/or fourth infection by age. Based on the relationship between FOI and % of people getting their second infection, we determine the FOI in Gampaha assuming that this region sees 77% of the number of second infections in Colombo. Based on the FOI we calculate the % of people experiencing their first or second infection each year, which is the recorded outcome in the clinical trial. using these estimated infection probabilities we can determine the trial sample size. 