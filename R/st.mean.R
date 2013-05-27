st.mean <-
function(df=10, sd=1, skew=1)
{
sqrt(df) * gamma((df-1)/2) * (skew-1/skew) * sd/(gamma(df/2) * sqrt(pi))
}
