# iomlifetR

An implementation of the IOMLIFET excel spreadsheets in R.

The key function is `impact()` which provides estimates of life-years gained in the future from a popoulation-wide reduction in PM~2.5~ related risk of mortality. 


## Additional things to consider

- Life table and impact function to include an exponential increase in risk. At the moment, the function doesn't cope particularly well with Australian data with an 85+ open ended age group. The mortality rate in this group is so low that you end up with a large (and unrealistic) number of people still alive at age 105. **-- DONE 9/1/17**

- Add a facility to alter birth rates

- Add a facility to accommodate migration. This  might be tricky because it requires assumptions about the past exposure of migrants and the application of cessation lag to these migrants.

- Make the impact function flexible enough to  accommodate abridged life tables.

- Adding a function to calculates the present value of life years gained from the data frame output by impact(), a VSLY, a discount rate and a time horizon of interest. 
