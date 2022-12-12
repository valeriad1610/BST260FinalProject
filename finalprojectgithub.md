Final Project Code
================
Valeria Duran
2022-12-12

## Introduction

This dataset, ‘Hatred on Twitter During MeToo Movement’, pulled 807,174
MeToo-related tweets from September 28, 2018 to February 16, 2019. The
MeToo Movement calls attention to the frequency at which, primarily
women and girls, experience sexual assault and harassment. This dataset
contains information on ten indicators: a unique ID for each tweet,
tweet text, timestamp of tweet, like count, retweet count, location
mentioned by user wile tweeting, follower count, following count, total
user tweet count, and whether tweet belongs to hatred or non-hated
(binary variable). The spread of hateful messaging on twitter is
harmful, especially in the context of the MeToo Movement, which is meant
to amplify the voices and experiences of survivors of sexual assault.
Hateful messaging on twitter diluted the MeToo campaign (Lingdren,
2019). Thus, it is important to understand what types of accounts and
what types of tweets sentiments spread hateful messaging. This project
explores the question: how do account attributes (i.e., following count,
follower count) and sentiment of tweets predict the hatefulness (binary
yes/no) of a user.

## Exploratory Analysis

As seen in **Figures** **1** and **2**, the distributions of follower
and following count were heavily skewed to the right, thus outliers were
removed. To further explore account attributes, the correlation between
follower count and following count per user was calculated; as seen in
**Fig. 3**, these two attributes demonstrate a slight positive
correlation, with r= 0.57 and p\>2.2e-16. Because there is no evident
strong correlation between following and follower counts, we will
proceed with exploring follower and following counts. As seen in
**Figures** **4** and **5**, there is a scattered trend between account
attributes and proportion of hateful tweets per user. However, there are
many users which tweet either entirely hateful or entirely non-hateful
tweets (proportion = 1, 0, respectively). There is a sizable
concentration of users who tweeted a proportion 0.1-0.3 hateful tweets
and another sizable concentration of users that have never tweeted
hatefully. Thus, we proceed by dichotomizing hateful accounts; a user is
hateful if they have ever tweeted hatefully (proportion of hateful
tweets\>0). Additionally, average sentiment score tends to be lower in
hateful users, as seen in **Fig. 6**. Thus, we will proceed with models
which explore the impact of average sentiment, follower count, an
following count on account hatefulness towards the MeToo Movement.

## Methods

Each user is identified and grouped by their unique combination of
following count, followers count, and location. Sentiment of tweets are
derived from “afinn” in the tidtytext package. An average sentiment
score is assigned to each user (calculated by averaging the sentiment of
words per tweet per user).

We will explore the prediction of account hatefulness by average
sentiment of users and account attributes using the caret package train
function with logistic regression models and 80/20 ratio of train/test,
based on a random sample of 40,000 tweets. The full predictive logistic
regression model will contain 3 features: following counts, follower
counts, and average sentiment scores per user. We first create
predictive model 1 using only 1 feature: average sentiment of tweets per
user. Subsequently, in predictive model 2 we will add follower and
following counts, having 3 features in total. We will obtain the
accuracies predictive models 1 and 2 with confusion matrices. We will
compare the distributions of predictive probabilities generated from the
predictive models 1 and 2 using an overlayed histogram.

As secondary analysis, we will include only users which have tweeted
more than once on the MeToo movement. We conduct the two predictive
logistic regression models with the same features as previously
described in predictive models 1 and 2. We will obtain the accuracies of
each model with confusion matrices and compare distributions of
predictive probabilities generated from train models in an overlayed
histogram.

## Results

The dataset pulled from Kaggle contained one column, ‘text’, which
contained the entire text of a tweet, including links, stop words, NAs,
and emojis. String processing and the tidyverse package were utilized to
remove stop words, links, emojis, and NA values from the data set. Then,
text mining methods in the tidytext package were utilized to expand the
dataset to include one row per word, assign sentiment values, and
average sentiment score grouped by user. A dichotomous variable for
hateful/non-hateful user was created. Then, variables for the average
sentiment scores per user and hateful user were joined with the original
dataset to include account attributes needed for analysis. This data
wrangling was necessary for our analysis, given our need for one data
frame which included an average sentiment score per twitter user, the
outcome variable (I.e., hateful/non-hateful user) and important
predictive attributes (I.e., following and follower counts). The data
before and after wrangling is demonstrated in Appendix: Data Wrangling

This project aims to understand how account attributes, specifically
following and follower count predict the hatefulness (hateful or not
hateful) of a twitter account toward the MeToo movement. To build the
predictive models, we first evaluated the association between average
sentiment of a tweet and the hatefulness of an account. We set a
prediction rule of 0.5; this is because a 0.5 is a standard cutoff. This
model generated an accuracy of 0.835, the confusion matrix of the
results is seen in **Fig. 7**. The model accurately predicted all
non-hateful accounts; however, it generated no prediction of hateful
accounts. In fact, all probabilities generated were below 0.5, as seen
in red in **Fig. 9**. Further, we added following and follower counts to
the logistic regression prediction model; using the same prediction rule
of 0.5, the full model generates the same accuracy and results, as seen
in **Fig. 8**. **Fig. 9** compares the probability distributions for the
full and partial models, which demonstrates that the full model (blue)
has a slightly higher maximum probability (0.31) than that of the model
which only used average sentiment as a predictor (0.28). Despite the
changes in predictive probabilities in the full model, the accuracy and
results of the full model remained unchanged, as seen in **Fig. 8**.

In the secondary analysis, we excluded users which tweeted only once,
this left 1153 users. This exclusion allows us to have a better estimate
of the hatefulness of a user, given the outcome variable (hatefulness of
user yes/no) is informed by at least two tweets.

Using this subset of the data, we then fit a logistic regression model
with average sentiment as the single predictor. With a 0.5 prediction
rule, this model generated a 73% accuracy. As shown in **Fig. 10**, this
model accurately predicted all non-hateful users and was unable to
predict any hateful users. When following and follower counts are added
to the model, the accuracy and results do not change, as seen in **Fig.
11**. However, the predictive probabilities differ between the partial
(red) and full models (blue), as seen in the comparison of their
distributions in **Fig. 12**. The full model generates slightly greater
predictive probabilities in comparison to the partial model. However,
these probabilities do not surpass the 0.5 probability decision
threshold.

## Conclusion

The aim of the project was to understand the predictive power of account
attributes (following and following counts) and average sentiment of
users on the hatefulness of a user towards the MeToo Movement. We
ascertained average sentiment of tweets for every unique user and used
this feature in a logistic regression model to predict the probability
that a user is hateful towards the MeToo Movement. We then added 2 more
features to the predictive logistic regression model: follower and
following counts. In the secondary analysis, the same models were fit to
data which included only user that have tweeted more than once about the
MeToo Movement. All models generated accuracies above 70%; however, no
models accurately predicted the hatefulness of a user. All models
generated predictive probabilities of hateful users below 0.4. Together,
results suggest that average sentiment, follower and following counts do
not strongly predict hatefulness of a user.

A limitation of the project was that we did not have username’s of users
and thus had to group users by following count, follower count, and
location; meaning, we may have incorrectly categorized the same users as
different users if any of the aforementioned variables changed over
time. This may have made our prediction model less accurate.

If I had more time, I would ask a different question of the data: how
does sentiment of tweets change over time? Using the data from the
secondary analysis (users who have tweeted more than once), I would have
compared the sentiment of tweets chronologically, to see how this varied
by user. This would allow us to understand how public opinion of the
MeToo movement changed or didn’t change over time.

## References

Lindgren, S. (2019). Movement Mobilization in the Age of Hashtag
Activism: Examining the Challenge of Noise, Hate, and Disengagement in
the \#MeToo Campaign. Policy & Internet, 11(4), 418–438.
<https://doi.org/10.1002/poi3.212>

## Appendix

``` r
library(tidytext)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(stringr)
library(textdata)
library(tidyverse)
```

    ## ── Attaching packages
    ## ───────────────────────────────────────
    ## tidyverse 1.3.2 ──

    ## ✔ ggplot2 3.4.0     ✔ readr   2.1.3
    ## ✔ tibble  3.1.8     ✔ purrr   0.3.5
    ## ✔ tidyr   1.2.1     ✔ forcats 0.5.2
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

``` r
library(caret)
```

    ## Loading required package: lattice
    ## 
    ## Attaching package: 'caret'
    ## 
    ## The following object is masked from 'package:purrr':
    ## 
    ##     lift

``` r
library(ggplot2)
library(readr)

data2 <- read_csv("~/Library/CloudStorage/OneDrive-HarvardUniversity/MeTooHate.csv")
```

    ## Rows: 807174 Columns: 10
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr  (2): text, location
    ## dbl  (7): status_id, favorite_count, retweet_count, followers_count, friends...
    ## dttm (1): created_at
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
#######need to change to make it able to knit
data2$index <- 1:nrow(data2) 
set.seed(11+16+2022)
dat<-data2[sample(nrow(data2), 40000), ]

# grouping tweets by the user 'newid'
dat2 <- dat |> unite(newid, location, followers_count, friends_count, remove=FALSE)
dat2$newid <- as.factor(dat2$newid)

#removing links and stop words 
links <- "https://t.co/[A-Za-z\\d]+|&amp;"
dat3 <- dat2 |> filter(!is.na(text)) |>
  mutate(text = str_replace_all(text, links, "")) |>
  unnest_tokens(word, text, token = "tweets")
```

    ## Using `to_lower = TRUE` with `token = 'tweets'` may not preserve URLs.

``` r
dat3 <- dat3 |> filter(!word %in% stop_words$word) 
```

``` r
#getting average sentiment 
afinn<-get_sentiments("afinn")
summary<- dat3 |> inner_join(afinn, by = "word") |> 
  group_by(newid) |>
  summarize(avgsent=mean(value))

#joining back together for a table that has account characteristics and the average sentiment score grouped by user 
inner_join(summary, dat2, by= 'newid') |> nrow() #27661
```

    ## [1] 27661

``` r
senttable <- inner_join(summary, dat2, by= 'newid')

#remove columns not needed 
senttable <- senttable |>
  subset(select = -c(status_id, created_at, favorite_count, retweet_count,statuses_count))
```

``` r
senttable <- senttable |>
  group_by(newid) |>
  mutate(prop = mean(category)) 
```

``` r
senttable <- senttable |> 
  group_by(newid) |>
  mutate(ntweet=n())
```

``` r
senttable |>
  group_by(newid) |>
  ggplot(aes(x=followers_count)) + 
  geom_histogram() +
  xlab('# of Followers') +
  labs(caption = 'Figure 1. The distribution of Followers by User') +
  theme(plot.caption=element_text(hjust = 0, face = "italic"))
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](finalprojectgithub_files/figure-gfm/expoloratory%20part%201-1.png)<!-- -->

``` r
senttable |>
 group_by(newid) |>
  ggplot(aes(x=friends_count)) + 
  geom_histogram() +
  xlab('# of Accounts Following') +
  labs(caption = 'Figure 2. The distribution of Accounts Following by User') +
  theme(plot.caption=element_text(hjust = 0, face = "italic"))
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](finalprojectgithub_files/figure-gfm/expoloratory%20part%201-2.png)<!-- -->

``` r
followeriqr<- IQR(senttable$followers_count)
summary(senttable$followers_count)
```

    ##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    ##        0       97      504    49762     2677 42719326

``` r
504 + (1.5*followeriqr)
```

    ## [1] 4374

``` r
sum(senttable$followers_count>=4374)
```

    ## [1] 5245

``` r
sub <- senttable |>
  filter(followers_count<=4374) 

followingiqr<- IQR(senttable$friends_count)
summary(senttable$friends_count)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##       0     150     508    7720    1709  521250

``` r
508 + (1.5*followingiqr)
```

    ## [1] 2846.5

``` r
sum(senttable$friends_count>=2846) 
```

    ## [1] 4829

``` r
sub <- sub |>
  filter(friends_count<=2846) 
```

``` r
cor.test(sub$followers_count, sub$friends_count) #0.5652156, p-value < 2.2e-16
```

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  sub$followers_count and sub$friends_count
    ## t = 98.276, df = 20574, p-value < 2.2e-16
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.5558446 0.5744430
    ## sample estimates:
    ##       cor 
    ## 0.5652156

``` r
#ADD SCATTERPLOT WITH LINE
sub |>
  group_by(newid) |>
  ggplot(aes(x=followers_count, y=friends_count)) + 
  geom_point()+
  geom_smooth(method = "lm") + 
  xlab("# of Followers") +
  ylab("# of Following") +
  labs(caption = 'Figure 3. Follower count versus following count per user') +
  theme(plot.caption=element_text(hjust = 0, face = "italic"))
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](finalprojectgithub_files/figure-gfm/exploratory%20part%202-1.png)<!-- -->

``` r
sub |> 
  group_by(newid) |>
  ggplot(aes(followers_count, prop, size=ntweet)) +
  geom_point() +
  ylab('Proportion of Hateful Tweets') +
  xlab('# of Followers') + 
  labs(caption = 'Figure 4. The proportion of hateful tweets versus number of followers by user') +
  theme(plot.caption=element_text(hjust = 0, face = "italic"))
```

![](finalprojectgithub_files/figure-gfm/exploratory%20part%203-1.png)<!-- -->

``` r
sub |> 
  group_by(newid) |>
  ggplot(aes(friends_count, prop, size=ntweet)) +
  geom_point() +
  ylab('Proportion of Hateful Tweets') +
  xlab('# of Accounts Following') + 
  labs(caption = 'Figure 5. The proportion of hateful tweets versus number of following count by user') +
  theme(plot.caption=element_text(hjust = 0, face = "italic"))
```

![](finalprojectgithub_files/figure-gfm/exploratory%20part%203-2.png)<!-- -->

``` r
sub <- sub |> 
  group_by(newid) |>
  mutate(prop = mean(category == 1)) |>
  mutate(hateful=case_when(
    prop > 0 ~ 1,
    prop == 0 ~ 0
  )) 

sub$hateful <- as.factor(sub$hateful)

sub |>
  ggplot(aes(x=hateful, y=avgsent)) + 
  geom_boxplot() +
  xlab('Hateful User') +
  ylab('Avergae Sentiment of Tweets per User') +
  scale_x_discrete(labels=c("1" = "Yes", "0" = "No"))+
  labs(caption = 'Figure 6. Average aentiment among Hateful and nonhateful users') +
  theme(plot.caption=element_text(hjust = 0, face = "italic"))
```

![](finalprojectgithub_files/figure-gfm/making%20dichotomozed%20hateful%20user%20variable-1.png)<!-- -->

``` r
logit1 <- glm(hateful ~ avgsent + followers_count, data=sub, family='binomial')
summary(logit1)
```

    ## 
    ## Call:
    ## glm(formula = hateful ~ avgsent + followers_count, family = "binomial", 
    ##     data = sub)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -0.9366  -0.6947  -0.5980  -0.4759   2.3634  
    ## 
    ## Coefficients:
    ##                   Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)     -1.533e+00  2.466e-02 -62.180  < 2e-16 ***
    ## avgsent         -1.922e-01  1.077e-02 -17.840  < 2e-16 ***
    ## followers_count -2.009e-04  2.553e-05  -7.867 3.63e-15 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 19765  on 20575  degrees of freedom
    ## Residual deviance: 19360  on 20573  degrees of freedom
    ## AIC: 19366
    ## 
    ## Number of Fisher Scoring iterations: 4

### Data Wrangling: Before & After

``` r
head(data2,1) 
```

    ## # A tibble: 1 × 11
    ##   status_id text     created_at          favor…¹ retwe…² locat…³ follo…⁴ frien…⁵
    ##       <dbl> <chr>    <dttm>                <dbl>   <dbl> <chr>     <dbl>   <dbl>
    ## 1   1.05e18 Entitle… 2018-09-30 01:17:15       5       1 McAlle…    2253    2303
    ## # … with 3 more variables: statuses_count <dbl>, category <dbl>, index <int>,
    ## #   and abbreviated variable names ¹​favorite_count, ²​retweet_count, ³​location,
    ## #   ⁴​followers_count, ⁵​friends_count

``` r
head(sub,1) 
```

    ## # A tibble: 1 × 11
    ## # Groups:   newid [1]
    ##   newid        avgsent text  locat…¹ follo…² frien…³ categ…⁴  index  prop ntweet
    ##   <fct>          <dbl> <chr> <chr>     <dbl>   <dbl>   <dbl>  <int> <dbl>  <int>
    ## 1 - GAME OVER…      -1 Moth… - GAME…       7      68       0 668030     0      1
    ## # … with 1 more variable: hateful <fct>, and abbreviated variable names
    ## #   ¹​location, ²​followers_count, ³​friends_count, ⁴​category

``` r
dsub<- sub[!duplicated(sub$newid), ]
y <- as.factor(dsub$hateful)
test_index <- createDataPartition(y, times = 1, p = 0.2, list = FALSE)

test_set <- dsub[test_index, ]
train_set <- dsub[-test_index, ]

nrow(train_set) # 14933
```

    ## [1] 14933

``` r
nrow(test_set) # 3735
```

    ## [1] 3735

``` r
test_set$hateful <- as.factor(test_set$hateful)
train_set$hateful <- as.factor(train_set$hateful)
lm_fit <- mutate(train_set, y = hateful == 1)   
lm_fit2<- glm(y ~ avgsent, data = lm_fit, family='binomial') 

#PREDICTION RULE OF 0.5
p_hat <- predict(lm_fit2, test_set, type='response')
y_hat <- ifelse(p_hat > 0.5, '1', '0') |> as.factor()
hist(p_hat)
```

![](finalprojectgithub_files/figure-gfm/using%20average%20user%20sentiment%20to%20predict%20hateful%20accounte-1.png)<!-- -->

``` r
max(p_hat)
```

    ## [1] 0.3133365

``` r
sd(p_hat)
```

    ## [1] 0.0484378

``` r
confusionMatrix(y_hat, test_set$hateful)$overall[["Accuracy"]] #0.8353414
```

    ## Warning in confusionMatrix.default(y_hat, test_set$hateful): Levels are not in
    ## the same order for reference and data. Refactoring data to match.

    ## [1] 0.8353414

``` r
cm1all <- confusionMatrix(data=y_hat, reference = test_set$hateful)
```

    ## Warning in confusionMatrix.default(data = y_hat, reference = test_set$hateful):
    ## Levels are not in the same order for reference and data. Refactoring data to
    ## match.

``` r
cm1all
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction    0    1
    ##          0 3120  615
    ##          1    0    0
    ##                                           
    ##                Accuracy : 0.8353          
    ##                  95% CI : (0.8231, 0.8471)
    ##     No Information Rate : 0.8353          
    ##     P-Value [Acc > NIR] : 0.5108          
    ##                                           
    ##                   Kappa : 0               
    ##                                           
    ##  Mcnemar's Test P-Value : <2e-16          
    ##                                           
    ##             Sensitivity : 1.0000          
    ##             Specificity : 0.0000          
    ##          Pos Pred Value : 0.8353          
    ##          Neg Pred Value :    NaN          
    ##              Prevalence : 0.8353          
    ##          Detection Rate : 0.8353          
    ##    Detection Prevalence : 1.0000          
    ##       Balanced Accuracy : 0.5000          
    ##                                           
    ##        'Positive' Class : 0               
    ## 

``` r
Reference <- factor(c(0, 0, 1, 1))
Prediction <- factor(c(0, 1, 0, 1))
y <- c(3120, 0, 615, 0)
df <- data.frame(Reference, Prediction, y)

ggplot(data =  df, mapping = aes(x = Reference, y = Prediction)) +
  geom_tile(aes(fill = y), colour = "white") +
  geom_text(aes(label = sprintf("%1.0f", y)), vjust = 1) +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_bw() + theme(legend.position = "none") +
  labs(caption = 'Figure 7. Actual hateful account vs. non-hateful account predictions using partial model') +
  theme(plot.caption=element_text(hjust = 0, face = "italic"))
```

![](finalprojectgithub_files/figure-gfm/confusion%20matrix%20plot%20for%20cm1%20all%20ntweet%20counts-1.png)<!-- -->

``` r
fit_glm <- glm(hateful ~ avgsent + followers_count+friends_count, data=train_set, family = "binomial")
p_hat_glm <- predict(fit_glm, test_set, type="response")
hist(p_hat_glm)
```

![](finalprojectgithub_files/figure-gfm/adding%20following%20and%20follower-1.png)<!-- -->

``` r
max(p_hat_glm)
```

    ## [1] 0.3271555

``` r
sd(p_hat_glm)
```

    ## [1] 0.05026637

``` r
y_hat_glm <- factor(ifelse(p_hat_glm > 0.5, '1', '0')) 
confusionMatrix(y_hat_glm, test_set$hateful)$overall["Accuracy"] #0.8353
```

    ## Warning in confusionMatrix.default(y_hat_glm, test_set$hateful): Levels are not
    ## in the same order for reference and data. Refactoring data to match.

    ##  Accuracy 
    ## 0.8353414

``` r
cm2all <- confusionMatrix(data=y_hat_glm, reference = test_set$hateful)
```

    ## Warning in confusionMatrix.default(data = y_hat_glm, reference =
    ## test_set$hateful): Levels are not in the same order for reference and data.
    ## Refactoring data to match.

``` r
cm2all
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction    0    1
    ##          0 3120  615
    ##          1    0    0
    ##                                           
    ##                Accuracy : 0.8353          
    ##                  95% CI : (0.8231, 0.8471)
    ##     No Information Rate : 0.8353          
    ##     P-Value [Acc > NIR] : 0.5108          
    ##                                           
    ##                   Kappa : 0               
    ##                                           
    ##  Mcnemar's Test P-Value : <2e-16          
    ##                                           
    ##             Sensitivity : 1.0000          
    ##             Specificity : 0.0000          
    ##          Pos Pred Value : 0.8353          
    ##          Neg Pred Value :    NaN          
    ##              Prevalence : 0.8353          
    ##          Detection Rate : 0.8353          
    ##    Detection Prevalence : 1.0000          
    ##       Balanced Accuracy : 0.5000          
    ##                                           
    ##        'Positive' Class : 0               
    ## 

``` r
Reference <- factor(c(0, 0, 1, 1))
Prediction <- factor(c(0, 1, 0, 1))
y <- c(3120, 0, 615, 0)
df <- data.frame(Reference, Prediction, y)

ggplot(data =  df, mapping = aes(x = Reference, y = Prediction)) +
  geom_tile(aes(fill = y), colour = "white") +
  geom_text(aes(label = sprintf("%1.0f", y)), vjust = 1) +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_bw() + theme(legend.position = "none") +
  labs(caption = 'Figure 8. Actual hateful account vs. non-hateful account predictions using full model') +
  theme(plot.caption=element_text(hjust = 0, face = "italic"))
```

![](finalprojectgithub_files/figure-gfm/confusion%20matrix%20plot%20for%20cm2%20all%20ntweet%20counts-1.png)<!-- -->

``` r
str(p_hat)
```

    ##  Named num [1:3735] 0.162 0.168 0.198 0.198 0.124 ...
    ##  - attr(*, "names")= chr [1:3735] "1" "2" "3" "4" ...

``` r
str(p_hat_glm)
```

    ##  Named num [1:3735] 0.17 0.151 0.178 0.209 0.13 ...
    ##  - attr(*, "names")= chr [1:3735] "1" "2" "3" "4" ...

``` r
dfp_hat <- data.frame(p_hat)
dfp_hat_glm <- data.frame(p_hat_glm)

dfp_hat_glm <- dfp_hat_glm |>
  rename("p" = "p_hat_glm")

dfp_hat <- dfp_hat |>
  rename("p" = "p_hat")

dfp <- rbind(dfp_hat, dfp_hat_glm)


ggplot(dfp, aes(x= p)) + 
  geom_histogram(data = dfp_hat, fill = "red", alpha = 0.2) + 
  geom_histogram(data = dfp_hat_glm, fill = "blue", alpha = 0.4) +
  theme() + 
   labs(caption = 'Figure 9. Comparison of probability distributions for hateful account predicitons') + 
    theme(plot.caption=element_text(hjust = 0, face = "italic")) + scale_color_manual(name='Histogram',
                     breaks=c('Partial Model', 'Full Model', 'Overlayed'),
                     values=c('Overalyed'='purple', 'Full model'='blue', 'Partial Model'= 'red')) +
                       theme(legend.title=element_text(size=20), legend.text=element_text(size=14))
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](finalprojectgithub_files/figure-gfm/how%20did%20probability%20change?-1.png)<!-- -->

``` r
########ask why legend isnt working 
```

``` r
table(dsub$ntweet>=2,dsub$hateful==1)
```

    ##        
    ##         FALSE  TRUE
    ##   FALSE 14757  2758
    ##   TRUE    840   313

``` r
sub2 <- sub |>
  filter(ntweet>=2)
nrow(sub2) 
```

    ## [1] 3061

``` r
dsub2<- sub2[!duplicated(sub2$newid), ]
nrow(dsub2)
```

    ## [1] 1153

``` r
y <- as.factor(dsub2$hateful)
test_index <- createDataPartition(y, times = 1, p = 0.2, list = FALSE)

test_set <- dsub2[test_index, ]
train_set <- dsub2[-test_index, ]

class(dsub2$hateful)
```

    ## [1] "factor"

``` r
test_set$hateful <- as.factor(test_set$hateful)
train_set$hateful <- as.factor(train_set$hateful)
lm_fit <- mutate(train_set, y = hateful == 1)   
lm_fit2<- lm(y ~ avgsent, data = lm_fit) 

#PREDICTION RULE OF 0.5
p_hat <- predict(lm_fit2, test_set)
hist(p_hat)
```

![](finalprojectgithub_files/figure-gfm/how%20does%20it%20change%20when%20we%20only%20take%20those%20that%20have%20tweeted%20more%20than%20once?-1.png)<!-- -->

``` r
y_hat <- ifelse(p_hat > 0.5, '1', '0') |> as.factor()
confusionMatrix(y_hat, test_set$hateful)$overall[["Accuracy"]] #0.7273
```

    ## Warning in confusionMatrix.default(y_hat, test_set$hateful): Levels are not in
    ## the same order for reference and data. Refactoring data to match.

    ## [1] 0.7272727

``` r
cm1 <- confusionMatrix(data=y_hat, reference = test_set$hateful)
```

    ## Warning in confusionMatrix.default(data = y_hat, reference = test_set$hateful):
    ## Levels are not in the same order for reference and data. Refactoring data to
    ## match.

``` r
cm1
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction   0   1
    ##          0 168  63
    ##          1   0   0
    ##                                          
    ##                Accuracy : 0.7273         
    ##                  95% CI : (0.665, 0.7836)
    ##     No Information Rate : 0.7273         
    ##     P-Value [Acc > NIR] : 0.5339         
    ##                                          
    ##                   Kappa : 0              
    ##                                          
    ##  Mcnemar's Test P-Value : 5.662e-15      
    ##                                          
    ##             Sensitivity : 1.0000         
    ##             Specificity : 0.0000         
    ##          Pos Pred Value : 0.7273         
    ##          Neg Pred Value :    NaN         
    ##              Prevalence : 0.7273         
    ##          Detection Rate : 0.7273         
    ##    Detection Prevalence : 1.0000         
    ##       Balanced Accuracy : 0.5000         
    ##                                          
    ##        'Positive' Class : 0              
    ## 

``` r
fit_glm <- glm(hateful ~ avgsent + followers_count+ friends_count, data=train_set, family = "binomial")
p_hat_glm <- predict(fit_glm, test_set, type="response")
hist(p_hat_glm)
```

![](finalprojectgithub_files/figure-gfm/how%20does%20it%20change%20when%20we%20only%20take%20those%20that%20have%20tweeted%20more%20than%20once?-2.png)<!-- -->

``` r
y_hat_glm <- factor(ifelse(p_hat_glm > 0.5, '1', '0')) 
confusionMatrix(y_hat_glm, test_set$hateful)$overall["Accuracy"] #0.7272727 
```

    ## Warning in confusionMatrix.default(y_hat_glm, test_set$hateful): Levels are not
    ## in the same order for reference and data. Refactoring data to match.

    ##  Accuracy 
    ## 0.7272727

``` r
cm2 <- confusionMatrix(data=y_hat_glm, reference = test_set$hateful)
```

    ## Warning in confusionMatrix.default(data = y_hat_glm, reference =
    ## test_set$hateful): Levels are not in the same order for reference and data.
    ## Refactoring data to match.

``` r
cm2
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction   0   1
    ##          0 168  63
    ##          1   0   0
    ##                                          
    ##                Accuracy : 0.7273         
    ##                  95% CI : (0.665, 0.7836)
    ##     No Information Rate : 0.7273         
    ##     P-Value [Acc > NIR] : 0.5339         
    ##                                          
    ##                   Kappa : 0              
    ##                                          
    ##  Mcnemar's Test P-Value : 5.662e-15      
    ##                                          
    ##             Sensitivity : 1.0000         
    ##             Specificity : 0.0000         
    ##          Pos Pred Value : 0.7273         
    ##          Neg Pred Value :    NaN         
    ##              Prevalence : 0.7273         
    ##          Detection Rate : 0.7273         
    ##    Detection Prevalence : 1.0000         
    ##       Balanced Accuracy : 0.5000         
    ##                                          
    ##        'Positive' Class : 0              
    ## 

``` r
cm1
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction   0   1
    ##          0 168  63
    ##          1   0   0
    ##                                          
    ##                Accuracy : 0.7273         
    ##                  95% CI : (0.665, 0.7836)
    ##     No Information Rate : 0.7273         
    ##     P-Value [Acc > NIR] : 0.5339         
    ##                                          
    ##                   Kappa : 0              
    ##                                          
    ##  Mcnemar's Test P-Value : 5.662e-15      
    ##                                          
    ##             Sensitivity : 1.0000         
    ##             Specificity : 0.0000         
    ##          Pos Pred Value : 0.7273         
    ##          Neg Pred Value :    NaN         
    ##              Prevalence : 0.7273         
    ##          Detection Rate : 0.7273         
    ##    Detection Prevalence : 1.0000         
    ##       Balanced Accuracy : 0.5000         
    ##                                          
    ##        'Positive' Class : 0              
    ## 

``` r
Reference <- factor(c(0, 0, 1, 1))
Prediction <- factor(c(0, 1, 0, 1))
y <- c(168, 0, 63, 0)
df <- data.frame(Reference, Prediction, y)

ggplot(data =  df, mapping = aes(x = Reference, y = Prediction)) +
  geom_tile(aes(fill = y), colour = "white") +
  geom_text(aes(label = sprintf("%1.0f", y)), vjust = 1) +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_bw() + theme(legend.position = "none") +
  labs(caption = 'Figure 10. Actual hateful account vs. non-hateful account predictions among users with more than one tweet using partial model') +
  theme(plot.caption=element_text(hjust = 0, face = "italic"))
```

![](finalprojectgithub_files/figure-gfm/confusion%20matrices-1.png)<!-- -->

``` r
cm2 
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction   0   1
    ##          0 168  63
    ##          1   0   0
    ##                                          
    ##                Accuracy : 0.7273         
    ##                  95% CI : (0.665, 0.7836)
    ##     No Information Rate : 0.7273         
    ##     P-Value [Acc > NIR] : 0.5339         
    ##                                          
    ##                   Kappa : 0              
    ##                                          
    ##  Mcnemar's Test P-Value : 5.662e-15      
    ##                                          
    ##             Sensitivity : 1.0000         
    ##             Specificity : 0.0000         
    ##          Pos Pred Value : 0.7273         
    ##          Neg Pred Value :    NaN         
    ##              Prevalence : 0.7273         
    ##          Detection Rate : 0.7273         
    ##    Detection Prevalence : 1.0000         
    ##       Balanced Accuracy : 0.5000         
    ##                                          
    ##        'Positive' Class : 0              
    ## 

``` r
ggplot(data =  df, mapping = aes(x = Reference, y = Prediction)) +
  geom_tile(aes(fill = y), colour = "white") +
  geom_text(aes(label = sprintf("%1.0f", y)), vjust = 1) +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_bw() + theme(legend.position = "none") +
  labs(caption = 'Figure 11. Actual hateful account vs. non-hateful account predictions among users with more than one tweet using full model') +
  theme(plot.caption=element_text(hjust = 0, face = "italic"))
```

![](finalprojectgithub_files/figure-gfm/confusion%20matrices-2.png)<!-- -->

``` r
str(p_hat)
```

    ##  Named num [1:231] 0.334 0.305 0.265 0.17 0.291 ...
    ##  - attr(*, "names")= chr [1:231] "1" "2" "3" "4" ...

``` r
str(p_hat_glm)
```

    ##  Named num [1:231] 0.368 0.187 0.131 0.204 0.125 ...
    ##  - attr(*, "names")= chr [1:231] "1" "2" "3" "4" ...

``` r
dfp_hat <- data.frame(p_hat)
dfp_hat_glm <- data.frame(p_hat_glm)

dfp_hat_glm <- dfp_hat_glm |>
  rename("p" = "p_hat_glm")

dfp_hat <- dfp_hat |>
  rename("p" = "p_hat")

dfp <- rbind(dfp_hat, dfp_hat_glm)


ggplot(dfp, aes(x= p)) + 
  geom_histogram(data = dfp_hat, fill = "red", alpha = 0.2) + 
  geom_histogram(data = dfp_hat_glm, fill = "blue", alpha = 0.4) +
  theme() + 
   labs(caption = 'Figure 12. Comparison of probability distributions for hateful account predicitons using accounts with more than 1 tweet') + 
    theme(plot.caption=element_text(hjust = 0, face = "italic")) + scale_color_manual(name='Histogram',
                     breaks=c('Partial Model', 'Full Model', 'Overlayed'),
                     values=c('Overalyed'='purple', 'Full model'='blue', 'Partial Model'= 'red')) 
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](finalprojectgithub_files/figure-gfm/how%20did%20probability%20change?%20part%202-1.png)<!-- -->
