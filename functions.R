get_accuracy <- function (data, model, label) {
  
  results <- matrix (NA, nrow = 50, ncol = 3)
  
  correct_bootstrap <- rep (NA, 10)
  
  for (i in 1:50) {
    
    lower <- i / 100
    
    upper <- 1 - lower
    
    predicted <- predict (model, data, type = 'prob')[,2]
    
    observed <- recode (pull (data, label), "yes" = 1, "no" = 0)
    
    actual <- tibble (
      predicted = predicted,
      observed = observed
    )
    
    correct_predictions <- actual %>%
      filter (predicted < lower | predicted > upper) %>% 
      filter ((predicted < lower & observed == 0) |
        (predicted > upper & observed == 1)) %>% 
      nrow ()
    
    predictions <- actual %>% 
       filter (predicted < lower | predicted > upper) %>% 
       nrow ()
    
    results [i] <- correct_predictions / predictions
  
    for (j in 1:1000) {

      data_bootstrap <- sample_n (data, size = nrow (data), replace = TRUE)

      predicted_bootstrap <- predict (model, data_bootstrap, type = 'prob')[,2]

      observed_bootstrap <- recode (pull (data_bootstrap, label), "yes" = 1, "no" = 0)

      bootstrap <- tibble (
        predicted = predicted_bootstrap,
        observed = observed_bootstrap
      )

      correct_predictions <- bootstrap %>%
        filter (predicted < lower | predicted > upper) %>%
        filter ((predicted < lower & observed == 0) |
          (predicted > upper & observed == 1)) %>%
        nrow ()

      predictions <- bootstrap %>%
        filter (predicted < lower | predicted > upper) %>%
        nrow ()

      correct_bootstrap [j] <- correct_predictions / predictions

    }

    results [i, 2] <- quantile (correct_bootstrap, probs = 0.025)
    results [i, 3] <- quantile (correct_bootstrap, probs = 0.975)

  }

  return (results)  
  
}

get_accuracy_ers_ats <- function (data) {
  
  correct_predictions <- data %>%
    filter (fvc_z_score > -1.645 & tlc_z_score > -1.645) %>% 
    nrow ()
    
  predictions <- data %>% 
    filter (fvc_z_score > -1.645) %>% 
    nrow ()
  
  accuracy <- correction_predictions / predictions
  
  return (accuracy)
  
}

get_cutoff <- function (model, data) {

  predicted <- predict (model, data, type = 'prob')[,2]
  
  observed <- recode (pull (data, restrictive), "yes" = 1, "no" = 0)
  
  suppressMessages (cutoff <- cutpointr (x = predicted, class = observed))
  
  return (cutoff$optimal_cutpoint)
  
}

get_false_positive <- function (model, data, label, threshold) {
  
  predicted <- predict (model, data, type = 'prob')[,2]
  observed <- recode (pull (data, label), "yes" = 1, "no" = 0)
  
  actual <- tibble (
    data = data,
    predicted = predicted,
    observed = observed
  )
    
  false_positive <- actual %>% 
    filter (predicted >= threshold & observed == 0) %>% 
    nrow () %>% 
    '/' (nrow (data))
    
  return (false_positive)    

}

get_false_positive_ers_ats <- function (data) {
  
  data <- data %>%  
    mutate (predicted = case_when (
      fvc_z_score < -1.645 ~ 1,
      TRUE ~ 0)
    ) %>% 
    mutate (observed = case_when (
      tlc_z_score < -1.645 ~ 1,
      TRUE ~ 0)
    )
  
  false_positive <- data %>% 
    filter (predicted == 1 & observed == 0) %>% 
    nrow () %>% 
    '/' (nrow (data))
  
  return (false_positive)
  
}

get_z_score <- function (age, height, sex, parameter, data) {
  
  mu <- get_mu (age, height, sex, parameter)
  sigma <- get_sigma (age, height, sex, parameter)
  nu <- get_nu (age, height, sex, parameter)
 
  z_score <- ((data / mu)^nu - 1) / (nu * sigma)
  
  return (z_score)
  
}

# GLI Global Reference Equations

get_mu <- function (age, height, sex, parameter) {
  
  if (parameter == "fev1") {
    
    if (sex == 1) {
      
      m_spline_lookup_table <- c (
        -0.076629683, -0.074514202, -0.074823324, -0.075695829, -0.072537805,  
        -0.061753529, -0.041712276, -0.01460646, 0.016886095, 0.049149418,
        0.078696905, 0.10346408, 0.123005703, 0.137322301, 0.146600993,
        0.151708117, 0.153621752, 0.153121563, 0.150830316, 0.147244378,
        0.142763215, 0.137710088, 0.132345912, 0.126834828, 0.121249596,
        0.11564601, 0.110070059, 0.104559617, 0.099136143, 0.093764967,
        0.088399368, 0.083000059, 0.077534064, 0.07197376, 0.06629601,
        0.060476196, 0.054484844, 0.048296381, 0.04188945, 0.035246319,
        0.028352392, 0.021195775, 0.01376715, 0.006068828, -0.001883729,
        -0.010075221, -0.018491388, -0.027118923, -0.035945394, -0.044959178,
        -0.054149393, -0.063505845, -0.073016305, -0.082655565, -0.092396128,
        -0.102212913, -0.112083045, -0.121985652, -0.131901696, -0.141813804,
        -0.151706132, -0.161564228, -0.171374918, -0.181126344, -0.190810569,
        -0.200422931, -0.209959313, -0.219416006, -0.228789679, -0.238077342,
        -0.247276324, -0.256384243, -0.265398986, -0.274318682, -0.283141688,
        -0.291866568, -0.300492075, -0.309017518, -0.317444317, -0.325774525,
        -0.334010132, -0.342153068, -0.350205204, -0.358168358, -0.366044293,
        -0.373834721, -0.381541304, -0.389165781, -0.396709555, -0.40417391
      )
      
      m_spline <- m_spline_lookup_table [age - 5]
      
      mu <- exp (-11.399108 + 2.462664 * log (height) - 0.011394 * log (age) + m_spline)
      
    } else {
      
      m_spline_lookup_table <- c (
        -0.138890127, -0.123387223, -0.106500168, -0.0867585, -0.063702899,
        -0.038030844, -0.011173899, 0.014961115, 0.039091821, 0.060315153,
        0.078122708, 0.092535584, 0.103844979, 0.112327738, 0.118398662,
        0.122576169, 0.125281242, 0.126849083, 0.127538078, 0.12755072,
        0.127051801, 0.126175989, 0.125010063, 0.123561972, 0.121829639,
        0.119813578, 0.117516272, 0.114940426, 0.112074018, 0.108897938,
        0.105397742, 0.10156304, 0.097386774, 0.092864634, 0.087995951,
        0.082786527, 0.077243858, 0.071375788, 0.065190389, 0.058695875,
        0.051900526, 0.044812637, 0.037444612, 0.029820512, 0.021964617, 
        0.013899103, 0.005644251, -0.002781368, -0.011360738, -0.020078246,
        -0.028919554, -0.037871298, -0.046917484, -0.056040141, -0.065222747,
        -0.074450222, -0.083708798, -0.0929859, -0.102270048, -0.111550754, 
        -0.120818444, -0.130064373, -0.139280562, -0.148459734, -0.157595268,
        -0.166681133, -0.175711838, -0.184682384, -0.193588228, -0.202425244,
        -0.211189692, -0.219878185, -0.228487664, -0.237015369, -0.245458821,
        -0.253815793, -0.262084298, -0.270263347, -0.278353995, -0.286357579,
        -0.294275416, -0.302108801, -0.309859007, -0.317527283, -0.325114858,
        -0.33262294, -0.340052714, -0.347405345, -0.354681976, -0.361883731
      )
      
      m_spline <- m_spline_lookup_table [age - 5]

      mu <- exp (-10.901689 + 2.385928 * log (height) - 0.076386 * log (age) + m_spline)
   
    }        
    
  }
  
  if (parameter == "fvc") {
    
    if (sex == 1) {
      
      m_spline_lookup_table <- c (
        -0.040465397, -0.042835353, -0.046580545, -0.04979453, -0.050163799,
        -0.046113283, -0.036958928, -0.02339395, -0.006386638, 0.012358413, 
        0.031070313, 0.048548818, 0.063991603, 0.076873763, 0.086919473,
        0.094355241, 0.099549029, 0.102808364, 0.104413805, 0.104692301, 
        0.103936201, 0.102385718, 0.100238459, 0.097639652, 0.094686887, 
        0.091460691, 0.088029015, 0.084449302, 0.080763297, 0.076969355, 
        0.073054326, 0.069007839, 0.064821851, 0.060490255, 0.056008528, 
        0.05137049, 0.046566513, 0.041588509, 0.036430171, 0.031086714, 
        0.025554642, 0.019831568, 0.013916223, 0.007815566, 0.001546046,
        -0.004876606, -0.011437984, -0.018124886, -0.024925196, -0.031827787,
        -0.038822425, -0.045899697, -0.053049639, -0.060256395, -0.067503357,
        -0.074775414, -0.082058811, -0.089341034, -0.096610696, -0.103857437, 
        -0.111071836, -0.118245329, -0.125370135, -0.132439298, -0.139448577,
        -0.146395987, -0.153279838, -0.160098642, -0.166851098, -0.173536074,
        -0.180152592, -0.186699817, -0.193177044, -0.199583686, -0.205919265,
        -0.212183402, -0.21837581, -0.224496486, -0.230546559, -0.236527485,
        -0.242440679, -0.248287512, -0.254069314, -0.259787376, -0.265442953,
        -0.271037262, -0.276571486, -0.282046886, -0.287464428, -0.292824984
      )
      
      m_spline <- m_spline_lookup_table [age - 5]

      mu <- exp (-12.629131 + 2.727421 * log (height) + 0.009174 * log (age) + m_spline)
      
    } else {
      
      m_spline_lookup_table <- c (
        -0.088232796, -0.078846087, -0.069656094, -0.059555991, -0.048029429,
        -0.035157502, -0.021361643, -0.007268505, 0.00661104, 0.019773006, 
        0.031886451, 0.042790016, 0.05242463, 0.060776064, 0.067920702,
        0.074001808, 0.07914435, 0.083454792, 0.087025974, 0.089938741, 
        0.092262928, 0.094059159, 0.095370919, 0.096211588, 0.0965914,
        0.096522124, 0.096016555, 0.095087507, 0.093740296, 0.091976457,
        0.089799516, 0.087214661, 0.084228359, 0.080848058, 0.077084193,
        0.072955247, 0.068480368, 0.063677496, 0.058563454, 0.053154032,
        0.047464061, 0.041507489, 0.035301158, 0.028871645, 0.022245107,
        0.015445279, 0.008493719, 0.001410025, -0.00578797, -0.013084028,
        -0.020463354, -0.02791235, -0.035416563, -0.042960981, -0.050531888,
        -0.058116805, -0.065704385, -0.073284305, -0.080847178, -0.088384471,
        -0.095888429, -0.103352005, -0.110768804, -0.118133438, -0.125442549,
        -0.132693467, -0.139883802, -0.147011417, -0.15407441, -0.161071094,
        -0.167999977, -0.174859746, -0.181649256, -0.18836751, -0.195013653,
        -0.201586955, -0.208086806, -0.214513177, -0.220867259, -0.227150393,
        -0.233363887, -0.239509012, -0.245587009, -0.251599087, -0.257546423,
        -0.263430165, -0.269251433, -0.275011318, -0.280710886, -0.286351175
      )
      
      m_spline <- m_spline_lookup_table [age - 5]

      mu <- exp (-12.055901 + 2.621579 * log (height) - 0.035975 * log (age) + m_spline)
      
    }
    
  }
  
  if (parameter == "fev1/fvc") {
    
    if (sex == 1) {
      
      m_spline_lookup_table <- c (
        -0.032117366, -0.031175746, -0.028435729, -0.024978945, -0.020173925, 
        -0.012645459, -0.00269177, 0.008311922, 0.018847292, 0.027821511,
        0.034654166, 0.039270883, 0.041939718, 0.043033726, 0.042939412, 
        0.041976643, 0.040417461, 0.038471711, 0.036290841, 0.033981428,
        0.031620531, 0.029262244, 0.026942503, 0.0246844, 0.022501808,
        0.020400952, 0.018377747, 0.016418939, 0.014508933, 0.012630404, 
        0.010775431, 0.008942585, 0.007133684, 0.005349458, 0.003591451,
        0.001860799, 0.000159386, -0.001511981, -0.003151699, -0.004759005,
        -0.006334722, -0.007883047, -0.00940954, -0.010917616, -0.012407048,
        -0.013877272, -0.015329917, -0.016766695, -0.018188443, -0.019594596,
        -0.020984153, -0.022356371, -0.023710902, -0.025047792, -0.026367909,
        -0.027671649, -0.028959947, -0.030233855, -0.031494384, -0.032431297,
        -0.033977614, -0.035201913, -0.036415157, -0.037616377, -0.03880421,
        -0.03997743, -0.04113479, -0.04227549, -0.043399198, -0.044505761,
        -0.045595232, -0.046667758, -0.047985281, -0.048763867, -0.049788563,
        -0.050798558, -0.051794487, -0.052776879, -0.053746213, -0.054702992,
        -0.05564775, -0.056580978, -0.057503132, -0.058414597, -0.059315692
      )
      
      m_spline <- m_spline_lookup_table [age - 5]

      mu <- exp (1.022608 - 0.218592 * log (height) - 0.027586 * log (age) + m_spline)
      
    } else {
      
      m_spline_lookup_table <- c (
        -0.044784354, -0.040463834, -0.03366171, -0.025586344, -0.01629574, 
        -0.005527463, 0.00580655, 0.016396842, 0.025252761, 0.031991009,
        0.036681939, 0.039611046, 0.041121214, 0.041523375, 0.041067832, 
        0.039979943, 0.038451661, 0.036636994, 0.034650922, 0.032567993,
        0.030436775, 0.028293127, 0.026159498, 0.024042924, 0.021945143, 
        0.019866836, 0.017808159, 0.015769248, 0.01375205, 0.01175805, 
        0.009795055, 0.007873685, 0.006004177, 0.004197303, 0.002464665,
        0.000816928, -0.000738675, -0.002199378, -0.00356516, -0.004838988,
        -0.006024997, -0.007129056, -0.008157648, -0.00911827, -0.010019727,
        -0.01087146, -0.011682432, -0.012459748, -0.013208761, -0.01393335,
        -0.014636918, -0.01532301, -0.015995898, -0.016659657, -0.01731823,
        -0.017975858, -0.018635432, -0.019298747, -0.019965877, -0.020635719, 
        -0.021307311, -0.021979595, -0.022651077, -0.023319829, -0.023984162,
        -0.024643034, -0.025295621, -0.025941356, -0.026579634, -0.027209901,
        -0.02783161, -0.028444489, -0.02904843, -0.029643611, -0.030230349,
        -0.030809071, -0.031380217, -0.031944124, -0.032501094, -0.033051396,
        -0.033595289, -0.034132946, -0.034664545, -0.035190276, -0.035710262
      )
      
      m_spline <- m_spline_lookup_table [age - 5]
      
      mu <- exp (0.9189568 - 0.1840671 * log (height) - 0.0461306 * log (age) + m_spline)
      
    }
    
  }
  
  if (parameter == "tlc") {
    
    if (sex == 1) {
      
      m_spline_lookup_table <- c (
        -0.084004394, -0.077960516, -0.071891043, -0.064588939, -0.055253542,
        -0.043587273, -0.029831587, -0.014709166, 0.000810758, 0.015785361,
        0.029433479, 0.041345825, 0.051346709, 0.059402186, 0.065638854,
        0.07019625, 0.073248057, 0.075002814, 0.075639347, 0.075323655,
        0.074229373, 0.072505502, 0.07027641, 0.067651297, 0.064726805,
        0.061584094, 0.058291245, 0.054906123, 0.051478833, 0.048051669,
        0.044660156, 0.041334056, 0.038092762, 0.034942205, 0.031885982,
        0.028926824, 0.026066728, 0.023306502, 0.020635958, 0.018036808,
        0.015492571, 0.012988674, 0.010512224, 0.008051817, 0.005596783,
        0.003135006, 0.000655001, -0.001853377, -0.004399081, -0.006990001,
        -0.009633081, -0.012334363, -0.015095595, -0.017912606, -0.020781005,
        -0.023696719, -0.026655957, -0.029655197, -0.03269116, -0.035760789,
        -0.038861137, -0.041987218, -0.045132265, -0.04828996, -0.051454492,
        -0.054620514, -0.057783109, -0.060937754, -0.064080289, -0.06720689,
        -0.070314055, -0.073399205, -0.076460982, -0.079498238, -0.082509929
      )
      
      m_spline <- m_spline_lookup_table [age - 5]
    
      mu <- exp (
        -10.5861 +
        0.1433 * log (age) +
        2.3155 * log (height) + 
        m_spline
      )      
      
    } else {
      
      m_spline_lookup_table <- c (
        -0.142986411, -0.119293793, -0.098957527, -0.081042267, -0.064516286,
        -0.048601064, -0.033160129, -0.018372153, -0.004709956, 0.007526767,
        0.018221651, 0.027348949, 0.0350004, 0.041306895, 0.046398371, 
        0.050446188, 0.05360938, 0.056019383, 0.057786724, 0.059004667,
        0.059751654, 0.060093439, 0.060081409, 0.059757591, 0.059158399,
        0.058315508, 0.057254369, 0.055993037, 0.054547241, 0.05293121,
        0.051157852, 0.049238768, 0.047183842, 0.045001942, 0.042701244,
        0.040289299, 0.037773093, 0.035158688, 0.032447289, 0.029637619,
        0.026728884, 0.023720742, 0.020613229, 0.017406699, 0.014101822,
        0.010701668, 0.007212254, 0.003639404, -1.15E-05, -0.003735307,
        -0.007527423, -0.011383424, -0.015299215, -0.019270414, -0.023288842,
        -0.027345163, -0.031430832, -0.035538037, -0.039659628, -0.04378907,
        -0.047920383, -0.052048096, -0.056167207, -0.060273167, -0.064362048,
        -0.068430383, -0.072475005, -0.076493028, -0.080481817, -0.084438976,
        -0.088362322, -0.092249873, -0.096099829, -0.099910561, -0.103680596
      )      
      
      m_spline <- m_spline_lookup_table [age - 5]
    
      mu <- exp (
        -10.1128 +
        0.1062 * log (age) +
        2.2259 * log (height) +
        m_spline
      ) 
      
    }
    
  }

  return (mu)
   
}

get_npv <- function (predicted, observed, cutoff) {
  
  data <- tibble (
    predicted = predicted,
    observed = observed
  )
  
  true_negatives <- data %>% 
    filter (predicted < cutoff & observed == 0) %>% 
    nrow ()
  
  false_negatives <- data %>% 
    filter (predicted < cutoff & observed == 1) %>% 
    nrow ()
  
  npv <- true_negatives / (true_negatives + false_negatives)
  
  return (npv)
  
}

get_npv_ers_ats <- function (data) {
  
  true_negatives <- data %>% 
    filter (fvc_z_score >= -1.645 & tlc_z_score >= -1.645) %>%
    nrow ()
    
  false_negatives <- data %>% 
    filter (fvc_z_score >= -1.645 & tlc_z_score < -1.645) %>%
    nrow ()      
  
  numerator <- true_negatives
  
  denominator <- true_negatives + false_negatives
  
  npv <- numerator / denominator
  
  lower <- exactci (numerator, denominator, 0.95)$conf.int [1]
  upper <- exactci (numerator, denominator, 0.95)$conf.int [2]

  return (
    list (
      npv = npv,
      lower = lower,
      upper = upper
    )
  )
  
}

# GLI Global Reference Equations

get_nu <- function (age, height, sex, parameter) {
  
  if (parameter == "fev1") {
    
    if (sex == 1) {
      
      nu <- 1.22703
      
    } else {
      
      nu <- 1.21388      
      
    }        
    
  }
  
  if (parameter == "fvc") {
    
    if (sex == 1) {
      
      nu <- 0.9346
      
    } else {
      
      nu <- 0.899
      
    }
    
  }
  
  if (parameter == "fev1/fvc") {
    
    if (sex == 1) {
      
      nu <- 3.8243 - 0.3328 * log (age)
      
    } else {
      
      nu <- 6.6490  - 0.9920 * log (age)
      
    }
    
  }
  
  if (parameter == "tlc") {
    
    if (sex == 1) {
      
      nu <- 0.9337
      
    } else {
      
      nu <- 0.4636
      
    }
    
  }
  
  return (nu)

}

get_performance <- function (model, data, cutoff) {
  
  bootstraps <- 10
  
  predicted <- predict (model, data, type = 'prob')[,2]
  
  observed <- recode (pull (data, restrictive), "yes" = 1, "no" = 0)

  sensitivity_bootstrap <- rep (NA, bootstraps)
  specificity_bootstrap <- rep (NA, bootstraps)
  npv_bootstrap <- rep (NA, bootstraps)
  ppv_bootstrap <- rep (NA, bootstraps)
  roc_bootstrap <- rep (NA, bootstraps)
  pr_bootstrap <- rep (NA, bootstraps)
  sbs_bootstrap <- rep (NA, bootstraps)
  ici_bootstrap <- rep (NA, bootstraps)
  
  for (i in 1:bootstraps) {
    
    data_bootstrap <- sample_n (data, size = nrow (data), replace = TRUE)
    
    predicted_bootstrap <- predict (model, data_bootstrap, type = 'prob')[,2]
    
    observed_bootstrap <- recode (pull (data_bootstrap, restrictive), "yes" = 1, "no" = 0)

    sensitivity_bootstrap [i] <- get_sensitivity (predicted_bootstrap, observed_bootstrap, cutoff)
    specificity_bootstrap [i] <- get_specificity (predicted_bootstrap, observed_bootstrap, cutoff)
    npv_bootstrap [i] <- get_npv (predicted_bootstrap, observed_bootstrap, cutoff)
    ppv_bootstrap [i] <- get_ppv (predicted_bootstrap, observed_bootstrap, cutoff)

    roc_bootstrap [i] <- AUC (predicted_bootstrap, observed_bootstrap)
    pr_bootstrap [i] <- PRAUC (predicted_bootstrap, observed_bootstrap)
    sbs_bootstrap [i] <- sbrier (predicted_bootstrap, observed_bootstrap)
    ici_bootstrap [i] <- suppressMessages (ici (predicted_bootstrap, observed_bootstrap))

  }
  
  results <- list (
    "sensitivity" = get_sensitivity (predicted, observed, cutoff),
    "sensitivity_025" = quantile (sensitivity_bootstrap, probs = 0.025, na.rm = TRUE)[[1]],
    "sensitivity_975" = quantile (sensitivity_bootstrap, probs = 0.975, na.rm = TRUE)[[1]],    
    "specificity" = get_specificity (predicted, observed, cutoff),
    "specificity_025" = quantile (specificity_bootstrap, probs = 0.025, na.rm = TRUE)[[1]],
    "specificity_975" = quantile (specificity_bootstrap, probs = 0.975, na.rm = TRUE)[[1]],    
    "npv" = get_npv (predicted, observed, cutoff),
    "npv_025" = quantile (npv_bootstrap, probs = 0.025, na.rm = TRUE)[[1]],
    "npv_975" = quantile (npv_bootstrap, probs = 0.975, na.rm = TRUE)[[1]],    
    "ppv" = get_ppv (predicted, observed, cutoff),
    "ppv_025" = quantile (ppv_bootstrap, probs = 0.025, na.rm = TRUE)[[1]],
    "ppv_975" = quantile (ppv_bootstrap, probs = 0.975, na.rm = TRUE)[[1]],    
    "roc" = AUC (predicted, observed),
    "roc_025" = quantile (roc_bootstrap, probs = 0.025, na.rm = TRUE)[[1]],
    "roc_975" = quantile (roc_bootstrap, probs = 0.975, na.rm = TRUE)[[1]],    
    "pr" = PRAUC (predicted, observed),
    "pr_025" = quantile (pr_bootstrap, probs = 0.025, na.rm = TRUE)[[1]],
    "pr_975" = quantile (pr_bootstrap, probs = 0.975, na.rm = TRUE)[[1]],    
    "ici" = suppressMessages (ici (predicted, observed)),
    "ici_025" = quantile (ici_bootstrap, probs = 0.025, na.rm = TRUE)[[1]],
    "ici_975" = quantile (ici_bootstrap, probs = 0.975, na.rm = TRUE)[[1]],       
    "sbs" = sbrier (predicted, observed),
    "sbs_025" = quantile (sbs_bootstrap, probs = 0.025, na.rm = TRUE)[[1]],
    "sbs_975" = quantile (sbs_bootstrap, probs = 0.975, na.rm = TRUE)[[1]]
  )
  
  return (results)
  
}

get_performance_ers_ats <- function (data) {
  
  bootstraps <- 10
  
  predicted <- data %>% 
    mutate (predicted = case_when (
      fvc_z_score < -1.645 ~ 1,
      fvc_z_score >= -1.645 ~ 0)
    ) %>% 
    pull (predicted)

  observed <- data %>% 
    mutate (observed = case_when (
      restrictive == "yes" ~ 1,
      restrictive == "no" ~ 0)
    ) %>% 
    pull (observed)
  
  sensitivity_bootstrap <- rep (NA, bootstraps)
  specificity_bootstrap <- rep (NA, bootstraps)
  npv_bootstrap <- rep (NA, bootstraps)
  ppv_bootstrap <- rep (NA, bootstraps)

  roc_bootstrap <- rep (NA, bootstraps)
  pr_bootstrap <- rep (NA, bootstraps)
  sbs_bootstrap <- rep (NA, bootstraps)

  for (i in 1:bootstraps) {
    
    data_bootstrap <- sample_n (data, size = nrow (data), replace = TRUE)
  
    predicted_bootstrap <- data_bootstrap %>% 
      mutate (predicted = case_when (
        fvc_z_score < -1.645 ~ 1,
        fvc_z_score >= -1.645 ~ 0)
      ) %>% 
      pull (predicted)

    observed_bootstrap <- data_bootstrap %>% 
      mutate (observed = case_when (
        restrictive == "yes" ~ 1,
        restrictive == "no" ~ 0)
      ) %>% 
      pull (observed)
    
    sensitivity_bootstrap [i] <- get_sensitivity_ers_ats (data_bootstrap)[[1]]
    specificity_bootstrap [i] <- get_specificity_ers_ats (data_bootstrap)[[1]]
    npv_bootstrap [i] <- get_npv_ers_ats (data_bootstrap)[[1]]
    ppv_bootstrap [i] <- get_ppv_ers_ats (data_bootstrap) [[1]]

    roc_bootstrap [i] <- AUC (predicted_bootstrap, observed_bootstrap)
    pr_bootstrap [i] <- PRAUC (predicted_bootstrap, observed_bootstrap)
    sbs_bootstrap [i] <- sbrier (predicted_bootstrap, observed_bootstrap)

  }
  
  results <- list (
    "sensitivity" = get_sensitivity_ers_ats (data)[[1]],
    "sensitivity_025" = quantile (sensitivity_bootstrap, probs = 0.025, na.rm = TRUE)[[1]],
    "sensitivity_975" = quantile (sensitivity_bootstrap, probs = 0.975, na.rm = TRUE)[[1]],    
    "specificity" = get_specificity_ers_ats (data)[[1]],
    "specificity_025" = quantile (specificity_bootstrap, probs = 0.025, na.rm = TRUE)[[1]],
    "specificity_975" = quantile (specificity_bootstrap, probs = 0.975, na.rm = TRUE)[[1]],    
    "npv" = get_npv_ers_ats (data)[[1]],
    "npv_025" = quantile (npv_bootstrap, probs = 0.025, na.rm = TRUE)[[1]],
    "npv_975" = quantile (npv_bootstrap, probs = 0.975, na.rm = TRUE)[[1]],    
    "ppv" = get_ppv_ers_ats (data)[[1]],
    "ppv_025" = quantile (ppv_bootstrap, probs = 0.025, na.rm = TRUE)[[1]],
    "ppv_975" = quantile (ppv_bootstrap, probs = 0.975, na.rm = TRUE)[[1]],
    "roc" = AUC (predicted, observed),
    "roc_025" = quantile (roc_bootstrap, probs = 0.025, na.rm = TRUE)[[1]],
    "roc_975" = quantile (roc_bootstrap, probs = 0.975, na.rm = TRUE)[[1]],    
    "pr" = PRAUC (predicted, observed),
    "pr_025" = quantile (pr_bootstrap, probs = 0.025, na.rm = TRUE)[[1]],
    "pr_975" = quantile (pr_bootstrap, probs = 0.975, na.rm = TRUE)[[1]],    
    "sbs" = sbrier (predicted, observed),
    "sbs_025" = quantile (sbs_bootstrap, probs = 0.025, na.rm = TRUE)[[1]],
    "sbs_975" = quantile (sbs_bootstrap, probs = 0.975, na.rm = TRUE)[[1]]
  )
  
  return (results)
  
}

get_ppv <- function (predicted, observed, cutoff) {
  
  data <- tibble (
    predicted = predicted,
    observed = observed
  )
  
  true_positives <- data %>% 
    filter (predicted > cutoff & observed == 1) %>% 
    nrow ()
  
  false_positives <- data %>% 
    filter (predicted > cutoff & observed == 0) %>% 
    nrow ()
  
  ppv <- true_positives / (true_positives + false_positives)
  
  return (ppv)
  
}

get_ppv_ers_ats <- function (data) {
  
  true_positives <- data %>% 
    filter (fvc_z_score < -1.645 & tlc_z_score < -1.645) %>%
    nrow ()
    
  false_positives <- data %>% 
    filter (fvc_z_score < -1.645 & tlc_z_score >= -1.645) %>%
    nrow ()  
  
  numerator <- true_positives
  
  denominator <- true_positives + false_positives

  ppv <- numerator / denominator
  
  return (ppv)
  
}

get_precision <- function (predicted, observed, cutoff) {
  
  tp <- 0
  fp <- 0
  
  for (i in 1:length (observed)) {
    
    if (observed [i] == 1 & predicted [i] >= cutoff) tp <- tp + 1
    if (observed [i] == 0 & predicted [i] >= cutoff) fp <- fp + 1
    
  }
  
  precision <- tp / (tp + fp)
  
  if (tp == 0 & fp == 0) precision <- 1
  
  return (precision)
  
}

get_precision_lower <- function (predicted, observed, cutoff) {
  
  bootstraps <- 10
  
  precision_bootstrap <- rep (NA, bootstraps)
  
  data <- tibble (
    observed = observed,
    predicted = predicted
  )
  
  for (i in 1:bootstraps) {
    
    data_bootstrap <- sample_n (data, size = nrow (data), replace = TRUE)
    
    tp <- nrow (filter (data_bootstrap, observed == 1 & predicted >= cutoff))
    
    fp <- nrow (filter (data_bootstrap, observed == 0 & predicted >= cutoff))    
    
    if (tp == 0 & fp == 0) {
      
      precision_bootstrap [i] <- 1
      
    } else {
    
      precision_bootstrap [i] <- tp / (tp + fp)
      
    }
    
  }
  
  precision_lower <- quantile (precision_bootstrap, probs = 0.025)

  return (precision_lower)    
  
}

get_precision_upper <- function (predicted, observed, cutoff) {
  
  bootstraps <- 10
  
  precision_bootstrap <- rep (NA, bootstraps)
  
  data <- tibble (
    observed = observed,
    predicted = predicted
  )
  
  for (i in 1:bootstraps) {
    
    data_bootstrap <- sample_n (data, size = nrow (data), replace = TRUE)
    
    tp <- nrow (filter (data_bootstrap, observed == 1 & predicted >= cutoff))
    
    fp <- nrow (filter (data_bootstrap, observed == 0 & predicted >= cutoff))    
    
    if (tp == 0 & fp == 0) {
      
      precision_bootstrap [i] <- 1
      
    } else {
    
      precision_bootstrap [i] <- tp / (tp + fp)
      
    }    
  }
  
  precision_upper <- quantile (precision_bootstrap, probs = 0.975)

  return (precision_upper)    
  
}

get_equity <- function (data, model, cutoff, measure) {
  
  performance <- rep (NA, 5)
  
  for (i in 3:4) {
    
    predicted <- predict (model, filter (data, race == i), type = 'prob')[, 2]
    observed <- recode (pull (filter (data, race == i), 'restrictive'), "yes" = 1, "no" = 0)
    
    if (measure == "sensitivity")
      performance [i] <- get_sensitivity (predicted, observed, cutoff)
    
    if (measure == "specificity")
      performance [i] <- get_specificity (predicted, observed, cutoff)
    
    if (measure == "npv")
      performance [i] <- get_npv (predicted, observed, cutoff)
    
    if (measure == "ppv")
      performance [i] <- get_ppv (predicted, observed, cutoff)
    
    if (measure == "roc")
      performance [i] <-  AUC (predicted, observed)
      
    if (measure == "pr")
      performance [i] <- PRAUC (predicted, observed)
      
    if (measure == "sbs")
      performance [i] <- sbrier (predicted, observed)
    
  }
  
  ratio <- min (performance, na.rm = TRUE) / max (performance, na.rm = TRUE) 

  return (ratio)     
  
}

get_equity_lower <- function (data, model, cutoff, measure) {
  
  bootstraps <- 10
  
  ratio_bootstrap <- rep (NA, bootstraps)
  
  for (i in 1:bootstraps) {

    performance <- rep (NA, 5)    
        
    data_bootstrap <- sample_n (data, size = nrow (data), replace = TRUE)
    
    for (j in 3:4) {
    
      predicted <- predict (model, filter (data_bootstrap, race == j), type = 'prob')[, 2]
      observed <- recode (pull (filter (data_bootstrap, race == j), 'restrictive'), "yes" = 1, "no" = 0)
        
      if (measure == "sensitivity")
        performance [j] <- get_sensitivity (predicted, observed, cutoff)
      
      if (measure == "specificity")
        performance [j] <- get_specificity (predicted, observed, cutoff)      
      
      if (measure == "npv")
        performance [j] <- get_npv (predicted, observed, cutoff)
      
      if (measure == "ppv")
        performance [j] <- get_ppv (predicted, observed, cutoff)
      
      if (measure == "roc")
        performance [j] <-  AUC (predicted, observed)
      
      if (measure == "pr")
        performance [j] <- PRAUC (predicted, observed)
      
      if (measure == "sbs")
        performance [j] <- sbrier (predicted, observed)
        
    }
    
    ratio_bootstrap [i] <- min (performance, na.rm = TRUE) / max (performance, na.rm = TRUE)

  }
  
  lower <- quantile (ratio_bootstrap, probs = 0.025)[[1]]
  
  return (lower)
  
}

get_equity_upper <- function (data, model, cutoff, measure) {
  
  bootstraps <- 10
  
  ratio_bootstrap <- rep (NA, bootstraps)
  
  for (i in 1:bootstraps) {

    performance <- rep (NA, 5)    
        
    data_bootstrap <- sample_n (data, size = nrow (data), replace = TRUE)
    
    for (j in 3:4) {
    
      predicted <- predict (model, filter (data_bootstrap, race == j), type = 'prob')[, 2]
      observed <- recode (pull (filter (data_bootstrap, race == j), 'restrictive'), "yes" = 1, "no" = 0)
        
      if (measure == "sensitivity")
        performance [j] <- get_sensitivity (predicted, observed, cutoff)
      
      if (measure == "specificity")
        performance [j] <- get_specificity (predicted, observed, cutoff)      
      
      if (measure == "npv")
        performance [j] <- get_npv (predicted, observed, cutoff)
      
      if (measure == "ppv")
        performance [j] <- get_ppv (predicted, observed, cutoff)
      
      if (measure == "roc")
        performance [j] <-  AUC (predicted, observed)
      
      if (measure == "pr")
        performance [j] <- PRAUC (predicted, observed)
      
      if (measure == "sbs")
        performance [j] <- sbrier (predicted, observed)
        
    }
    
    ratio_bootstrap [i] <- min (performance, na.rm = TRUE) / max (performance, na.rm = TRUE)

  }
  
  upper <- quantile (ratio_bootstrap, probs = 0.975)[[1]]
  
  return (upper)
  
}

get_equity_ers_ats <- function (data, measure) {
  
  performance <- rep (NA, 5)
  
  for (i in 3:4) {
    
    if (measure == "sensitivity")
      performance [i] <- get_sensitivity_ers_ats (filter (data, race == i))
      
    if (measure == "specificity") 
      performance [i] <- get_specificity_ers_ats (filter (data, race == i))
      
    if (measure == "npv")
      performance [i] <- get_npv_ers_ats (filter (data, race == i))[[1]]
      
    if (measure == "ppv")       
      performance [i] <- get_ppv_ers_ats (filter (data, race == i))
    
    if (measure == "roc" | measure == "pr" | measure == "sbs") {
      
      predicted <- data %>%
        filter (race == i) %>% 
        mutate (predicted = case_when (
          fvc_z_score < -1.645 ~ 1,
          TRUE ~ 0)
        ) %>% 
        pull (predicted)
      
      observed <- data %>%
        filter (race == i) %>% 
        mutate (observed = case_when (
          tlc_z_score < -1.645 ~ 1,
          TRUE ~ 0)
        ) %>% 
        pull (observed)
      
      if (measure == "roc")
        performance [i] <- AUC (predicted, observed)
      
      if (measure == "pr")
        performance [i] <- PRAUC (predicted, observed)
      
      if (measure == "sbs")
        performance [i] <- sbrier (predicted, observed)
     
    }
    
  }
  
  ratio <- min (performance, na.rm = TRUE) / max (performance, na.rm = TRUE)

  return (ratio)   

}

get_equity_ers_ats_lower <- function (data, measure) {
  
  bootstraps <- 10
  
  ratio_bootstrap <- rep (NA, bootstraps)
  
  for (i in 1:bootstraps) {

    performance <- rep (NA, 5)    
        
    data_bootstrap <- sample_n (data, size = nrow (data), replace = TRUE)
    
    for (j in 3:4) {
    
      if (measure == "sensitivity")
        performance [j] <- get_sensitivity_ers_ats (filter (data_bootstrap, race == j))
      
      if (measure == "specificity") 
        performance [j] <- get_specificity_ers_ats (filter (data_bootstrap, race == j))
      
      if (measure == "npv")
        performance [j] <- get_npv_ers_ats (filter (data_bootstrap, race == j))[[1]]
      
      if (measure == "ppv")       
        performance [j] <- get_ppv_ers_ats (filter (data_bootstrap, race == j))
    
      if (measure == "roc" | measure == "pr" | measure == "sbs") {
      
        predicted <- data_bootstrap %>%
          filter (race == j) %>% 
          mutate (predicted = case_when (
            fvc_z_score < -1.645 ~ 1,
            TRUE ~ 0)
          ) %>% 
          pull (predicted)
      
        observed <- data_bootstrap %>%
          filter (race == j) %>% 
          mutate (observed = case_when (
            tlc_z_score < -1.645 ~ 1,
            TRUE ~ 0)
          ) %>% 
          pull (observed)
      
        if (measure == "roc")
          performance [j] <- AUC (predicted, observed)
      
        if (measure == "pr")
          performance [j] <- PRAUC (predicted, observed)
      
        if (measure == "sbs")
          performance [j] <- sbrier (predicted, observed)
        
      }
     
    }
    
    ratio_bootstrap [i] <- min (performance, na.rm = TRUE) / max (performance, na.rm = TRUE)
    
  }
  
  lower <- quantile (ratio_bootstrap, probs = 0.025)[[1]]
      
  return (lower)
  
}

get_equity_ers_ats_upper <- function (data, measure) {
  
  bootstraps <- 10
  
  ratio_bootstrap <- rep (NA, bootstraps)
  
  for (i in 1:bootstraps) {

    performance <- rep (NA, 5)    
        
    data_bootstrap <- sample_n (data, size = nrow (data), replace = TRUE)
    
    for (j in 3:4) {
    
      if (measure == "sensitivity")
        performance [j] <- get_sensitivity_ers_ats (filter (data_bootstrap, race == j))
      
      if (measure == "specificity") 
        performance [j] <- get_specificity_ers_ats (filter (data_bootstrap, race == j))
      
      if (measure == "npv")
        performance [j] <- get_npv_ers_ats (filter (data_bootstrap, race == j))[[1]]
      
      if (measure == "ppv")       
        performance [j] <- get_ppv_ers_ats (filter (data_bootstrap, race == j))
    
      if (measure == "roc" | measure == "pr" | measure == "sbs") {
      
        predicted <- data_bootstrap %>%
          filter (race == j) %>% 
          mutate (predicted = case_when (
            fvc_z_score < -1.645 ~ 1,
            TRUE ~ 0)
          ) %>% 
          pull (predicted)
      
        observed <- data_bootstrap %>%
          filter (race == j) %>% 
          mutate (observed = case_when (
            tlc_z_score < -1.645 ~ 1,
            TRUE ~ 0)
          ) %>% 
          pull (observed)
      
        if (measure == "roc")
          performance [j] <- AUC (predicted, observed)
      
        if (measure == "pr")
          performance [j] <- PRAUC (predicted, observed)
      
        if (measure == "sbs")
          performance [j] <- sbrier (predicted, observed)
        
      }
     
    }
    
    ratio_bootstrap [i] <- min (performance, na.rm = TRUE) / max (performance, na.rm = TRUE)
    
  }
  
  upper <- quantile (ratio_bootstrap, probs = 0.975)[[1]]
      
  return (upper)
  
}

get_sbs <- function (model, data) {
  
  predicted <- predict (model, data, type = 'prob')[,2]
  observed <- recode (pull (data, "restrictive"), "yes" = 1, "no" = 0)  
  
  sbs <- sbrier (predicted, observed)
  
  return (sbs)
  
}

get_sbs_lower <- function (model, data) {
  
  bootstraps <- 10
  
  sbs_bootstrap <- rep (NA, bootstraps)
  
  for (i in 1:bootstraps) {
    
    data_bootstrap <- sample_n (data, size = nrow (data), replace = TRUE)

    predicted <- predict (model, data_bootstrap, type = 'prob')[,2]
    observed <- recode (pull (data_bootstrap, "restrictive"), "yes" = 1, "no" = 0)  
  
    sbs_bootstrap [i] <- sbrier (predicted, observed)
    
  }
  
  sbs_lower <- quantile (sbs_bootstrap, probs = 0.025)[[1]]
  
  return (sbs_lower)
  
}

get_sbs_upper <- function (model, data) {
  
  bootstraps <- 10
  
  sbs_bootstrap <- rep (NA, bootstraps)
  
  for (i in 1:bootstraps) {
    
    data_bootstrap <- sample_n (data, size = nrow (data), replace = TRUE)

    predicted <- predict (model, data_bootstrap, type = 'prob')[,2]
    observed <- recode (pull (data_bootstrap, "restrictive"), "yes" = 1, "no" = 0)  
  
    sbs_bootstrap [i] <- sbrier (predicted, observed)
    
  }
  
  sbs_lower <- quantile (sbs_bootstrap, probs = 0.975)[[1]]
  
  return (sbs_lower)
  
}

get_sbs_ers_ats <- function (data) {
  
  predicted <- data %>% 
    mutate (predicted = case_when (
      fvc_z_score > -1.645 ~ 0,
      TRUE ~ 1)
    ) %>% 
    pull (predicted)
  
  observed <- data %>% 
    mutate (observed = case_when (
      tlc_z_score < -1.645 ~ 1,
      TRUE ~ 0)
    ) %>% 
    pull (observed)
  
  sbs <- sbrier (predicted, observed)
  
  return (sbs)
  
}

get_sbs_ers_ats_lower <- function (data) {
  
  bootstraps <- 10
  
  sbs_bootstrap <- rep (NA, bootstraps)
  
  for (i in 1:bootstraps) {
    
    data_bootstrap <- sample_n (data, size = nrow (data), replace = TRUE)

    predicted <- data_bootstrap %>% 
      mutate (predicted = case_when (
        fvc_z_score > -1.645 ~ 0,
        TRUE ~ 1)
      ) %>% 
      pull (predicted)
  
    observed <- data_bootstrap %>% 
      mutate (observed = case_when (
        tlc_z_score < -1.645 ~ 1,
        TRUE ~ 0)
      ) %>% 
      pull (observed)
  
    sbs_bootstrap [i] <- sbrier (predicted, observed)
    
  }
  
  sbs_lower <- quantile (sbs_bootstrap, probs = 0.025)[[1]]
  
  return (sbs_lower)
  
}

get_sbs_ers_ats_upper <- function (data) {
  
  bootstraps <- 10
  
  sbs_bootstrap <- rep (NA, bootstraps)
  
  for (i in 1:bootstraps) {
    
    data_bootstrap <- sample_n (data, size = nrow (data), replace = TRUE)

    predicted <- data_bootstrap %>% 
      mutate (predicted = case_when (
        fvc_z_score > -1.645 ~ 0,
        TRUE ~ 1)
      ) %>% 
      pull (predicted)
  
    observed <- data_bootstrap %>% 
      mutate (observed = case_when (
        tlc_z_score < -1.645 ~ 1,
        TRUE ~ 0)
      ) %>% 
      pull (observed)
  
    sbs_bootstrap [i] <- sbrier (predicted, observed)
    
  }
  
  sbs_upper <- quantile (sbs_bootstrap, probs = 0.975)[[1]]
  
  return (sbs_upper)
  
}

get_sensitivity <- function (predicted, observed, cutoff) {
  
  tp <- 0 # True positives
  fn <- 0 # False negatives
  
  for (i in 1:length (observed)) {
    
    if (observed [i] == 1 & predicted [i] > cutoff) tp <- tp + 1
    
    if (observed [i] == 1 & predicted [i] <= cutoff) fn <- fn + 1
    
  }

  sensitivity <- tp / (tp + fn)
  
  return (sensitivity)
  
}

get_sensitivity_ers_ats <- function (data) {
  
  true_positives <- data %>% 
    filter (fvc_z_score < -1.645 & tlc_z_score < -1.645) %>%
      nrow ()
    
  false_negatives <- data %>% 
    filter (fvc_z_score >= -1.645 & tlc_z_score < -1.645) %>%
    nrow ()
  
  numerator <- true_positives
  
  denominator <- true_positives + false_negatives
  
  sensitivity <- numerator / denominator
  
  return (sensitivity)
  
}

get_sensitivity_lower <- function (predicted, observed, cutoff) {
  
  bootstraps <- 10
  
  sensitivity_bootstrap <- rep (NA, bootstraps)
  
  data <- tibble (
    observed = observed,
    predicted = predicted
  )
  
  for (i in 1:bootstraps) {
    
    data_bootstrap <- sample_n (data, size = nrow (data), replace = TRUE)
    
    tp <- nrow (filter (data_bootstrap, observed == 1 & predicted > cutoff))
    
    fn <- nrow (filter (data_bootstrap, observed == 1 & predicted <= cutoff))    
    
    sensitivity_bootstrap [i] <- tp / (tp + fn)
    
  }
  
  sensitivity_lower <- quantile (sensitivity_bootstrap, probs = 0.025)

  return (sensitivity_lower)
  
}

get_sensitivity_upper <- function (predicted, observed, cutoff) {
  
  bootstraps <- 10
  
  sensitivity_bootstrap <- rep (NA, bootstraps)
  
  data <- tibble (
    observed = observed,
    predicted = predicted
  )
  
  for (i in 1:bootstraps) {
    
    data_bootstrap <- sample_n (data, size = nrow (data), replace = TRUE)
    
    tp <- nrow (filter (data_bootstrap, observed == 1 & predicted > cutoff))
    
    fn <- nrow (filter (data_bootstrap, observed == 1 & predicted <= cutoff)) 
  
    sensitivity_bootstrap [i] <- tp / (tp + fn)
    
  }
  
  sensitivity_upper <- quantile (sensitivity_bootstrap, probs = 0.975)

  return (sensitivity_upper)
  
}

# GLI Global Reference Equations

get_sigma <- function (age, height, sex, parameter) {
  
  if (parameter == "fev1") {
    
    if (sex == 1) {
      
      s_spline_lookup_table <- c (
        0.072337153, 0.049397165, 0.030805517, 0.015607502, 0.003078378,
        -0.00735457, -0.01612753, -0.023555217, -0.029871683, -0.035233793,
        -0.039754188, -0.043523323, -0.046599967, -0.049034207, -0.05087416,
        -0.052170268, -0.052970381, -0.053317212, -0.053251042, -0.052816226,
        -0.052053473, -0.050998097, -0.049680914, -0.048130148, -0.046372107,
        -0.044429906, -0.04232383, -0.040071737, -0.037689471, -0.035191432,
        -0.032590514, -0.029898173, -0.027124612, -0.02427894, -0.021369316,
        -0.018403848, -0.015391107, -0.01233885, -0.009253983, -0.006142663,
        -0.003010387, 0.000137935, 0.003297898, 0.006464856, 0.009633744,
        0.012799988, 0.015959512, 0.019108692, 0.022244298, 0.02536346,
        0.02846363, 0.031542545, 0.034598212, 0.037628925, 0.040633199,
        0.043609727, 0.046557367, 0.049475123, 0.052362129, 0.05521764,
        0.05804102, 0.060831727, 0.063589309, 0.066313411, 0.069004114,
        0.071661853, 0.074287076, 0.076880233, 0.079441771, 0.081972135,
        0.08447177, 0.086941116, 0.089380608, 0.091790681, 0.094171762,
        0.096524275, 0.098848639, 0.101145296, 0.103414826, 0.105657843,
        0.107874938, 0.110066682, 0.11223363, 0.114376316, 0.116495259,
        0.118590961, 0.120663908, 0.122714586, 0.124743431, 0.126750856
      )
      
      s_spline <- s_spline_lookup_table [age - 5]
      
      sigma <- exp (-2.256278 + 0.080729 * log (age) + s_spline)

    } else {
      
      s_spline_lookup_table <- c (
        0.105023469, 0.066955942, 0.039364403, 0.019945763, 0.006619247, 
        -0.002596684, -0.009152836, -0.014111829, -0.018161267, -0.02171957,
        -0.025063881, -0.028344845, -0.031637911, -0.034992853, -0.038398682, 
        -0.041798173 ,-0.045146233, -0.048405921, -0.051512011, -0.054392748, 
        -0.056994145, -0.059276379, -0.061207648, -0.062757208, -0.063901867, 
        -0.064625772, -0.064918883, -0.064776551, -0.064207446, -0.063227952,
        -0.061854231, -0.060101998, -0.057986465, -0.05552232, -0.052727422, 
        -0.049631436, -0.046263962, -0.042651664, -0.038818604, -0.034786534, 
        -0.030575146, -0.026202299, -0.02168616, -0.017048938, -0.012311671,
        -0.007493286, -0.002610817, 0.002320395, 0.007286528, 0.012275129,
        0.017274972, 0.022275989, 0.027270092, 0.032250877, 0.037212604, 
        0.042150101, 0.047058712, 0.051934248, 0.056772944, 0.061571423,
        0.066326654, 0.07103593, 0.07569683, 0.080307548, 0.084867772,
        0.089377572, 0.093837067, 0.098246421, 0.102605840, 0.106915564, 
        0.111175867, 0.11538705, 0.119549439, 0.123663383, 0.127729251,
        0.131747427, 0.135718311, 0.139642476, 0.143520879, 0.147354507, 
        0.151144316, 0.154891231, 0.158596149, 0.162259939, 0.165883443, 
        0.169467478, 0.173012837, 0.176520288, 0.179990577, 0.18342443
      )
      
      s_spline <- s_spline_lookup_table [age - 5]
      
      sigma <- exp (-2.364047 + 0.129402 * log (age) + s_spline)

    }        
    
  }
  
  if (parameter == "fvc") {
    
    if (sex == 1) {
    
      s_spline_lookup_table <- c (
        0.057281204, 0.038624578, 0.023607021, 0.011424295, 0.001471157,
        -0.006723165, -0.013517123, -0.019178381, -0.023909971, -0.027865889,
        -0.031165881, -0.033901236, -0.036128414, -0.037892547, -0.039233567,
        -0.040183109, -0.040768573, -0.041015495, -0.04094946, -0.040601584,
        -0.040001133, -0.03917367, -0.038141656, -0.036926358, -0.035548354,
        -0.034025774, -0.032374526, -0.030608619, -0.028740637, -0.026782855,
        -0.024746482, -0.022641419, -0.020476429, -0.018259282, -0.015996891,
        -0.013695933, -0.011363211, -0.009004838, -0.006626239, -0.004232229,
        -0.00182709, 0.000585371, 0.003001762, 0.005418772, 0.007833085,
        0.010241727, 0.012642058, 0.015031736, 0.017408682, 0.019771057,
        0.022117233, 0.024445772, 0.02675542, 0.029045133, 0.031314018, 
        0.033561296, 0.035786295, 0.037988433, 0.040167216, 0.042322224, 
        0.044453106, 0.046559571, 0.048641387, 0.05069838, 0.05273063,
        0.054738424, 0.05672206, 0.058681838, 0.060618061, 0.062531032, 
        0.064421054, 0.066288432, 0.068133468, 0.069956463, 0.071757718,
        0.07353753, 0.075296195, 0.077034028, 0.078751454, 0.080448924,
        0.082126875, 0.083785727, 0.085425888, 0.087047756, 0.088651712,
        0.090238127, 0.091807362, 0.093359781, 0.0948957, 0.096415417
      )
      
      s_spline <- s_spline_lookup_table [age - 5]
      
      sigma <- exp (-2.195595 + 0.068466 * log (age) + s_spline)
      
    } else {
      
      s_spline_lookup_table <- c (
        0.091010919, 0.05893963, 0.036122984, 0.020503561, 0.010163041,
        0.003261599, -0.001542841, -0.00523079, -0.008435044, -0.011550642,
        -0.014829047, -0.018365664, -0.022169975, -0.026239075, -0.030504927,
        -0.034849674, -0.039181104, -0.043426446, -0.047495, -0.051296751,
        -0.0547639, -0.057846269, -0.060506549, -0.062715991, -0.064454526,
        -0.065709309, -0.066473193, -0.066744244, -0.066533134, -0.065857825,
        -0.064736088, -0.063185257, -0.061222168, -0.058863112, -0.056128272,
        -0.053052068, -0.049668916, -0.046009776, -0.042102551, -0.037972429,
        -0.033642183, -0.029132437, -0.024464041, -0.0196622, -0.014750744,
        -0.009751093, -0.004682517, 0.000437644, 0.005593781, 0.010771837,
        0.015959154, 0.021144425, 0.026319268, 0.031477755, 0.036614607,
        0.041725073, 0.04680488, 0.051850189, 0.056857557, 0.061823895, 
        0.066746443, 0.071622733, 0.076450568, 0.08122829, 0.085955514, 
        0.090632196, 0.095258352, 0.099834049, 0.104359404, 0.108834577,
        0.113259766, 0.117635204, 0.121961154, 0.126237907, 0.130465777,
        0.1346451, 0.138776234, 0.142859703, 0.146896415, 0.150887308,
        0.154833293, 0.158735256, 0.162594055, 0.166410525, 0.170185478,
        0.173919703, 0.177613968, 0.181269017, 0.184885577, 0.188464355
      )
      
      s_spline <- s_spline_lookup_table [age - 5]
      
      sigma <- exp (-2.310148 + 0.120428 * log (age) + s_spline)

    }
    
  }
  
  if (parameter == "fev1/fvc") {
    
    if (sex == 1) {
      
      s_spline_lookup_table <- c (
        0.04501273, 0.044271869, 0.041741105, 0.036751685, 0.030798272, 
        0.026562318, 0.024563815, 0.024204644, 0.024432178, 0.0244019,
        0.023533943, 0.021612804, 0.018671858, 0.014840877, 0.010335984,
        0.005396413, 0.000227949, -0.005004279, -0.010165307, -0.01514917,
        -0.019858591, -0.024189964, -0.028050214, -0.031361205, -0.034067995,
        -0.036131842, -0.037538767, -0.038290925, -0.038399416, -0.037872532,
        -0.036716206, -0.034939952, -0.032560508, -0.029599422, -0.026080176,
        -0.022028226, -0.017468703, -0.01242601, -0.006929132, -0.001009467,
        0.005300503, 0.011965483, 0.01894945, 0.026220461, 0.033753141,
        0.041522129, 0.0495015, 0.057665355, 0.065988728, 0.074447997,
        0.083021379, 0.091689499, 0.100433447, 0.109233818, 0.11807297,
        0.126935006, 0.135804097, 0.144665683, 0.153505505, 0.162310495, 
        0.171069048, 0.179771615, 0.188410433, 0.196979232, 0.205472949,
        0.213887511, 0.22222013, 0.230468683, 0.238631437, 0.246707212, 
        0.254695343, 0.262595384, 0.270406884, 0.278129776, 0.285764342, 
        0.29331088, 0.300770291, 0.308143798, 0.315432574, 0.322637725, 
        0.329760361, 0.336801784, 0.343763346, 0.35064663, 0.357453286
      )
      
      s_spline <- s_spline_lookup_table [age - 5]
      
      sigma <- exp (-2.882025 + 0.068889 * log (age) + s_spline)

    } else {
      
      s_spline_lookup_table <- c (
        0.051845823, 0.043555518, 0.03604702, 0.028926881, 0.022262836, 
        0.016386321, 0.011249676, 0.00668324, 0.00254365, -0.00124848,
        -0.004746278, -0.007968582, -0.010921835, -0.013622212, -0.016093567,
        -0.018351395, -0.020407135, -0.022265797, -0.023925738, -0.025387753,
        -0.026649986, -0.02770541, -0.028543983, -0.02915278, -0.029520715,
        -0.029641288, -0.029507811, -0.029114184, -0.028453694, -0.027520198,
        -0.026307581, -0.024811479, -0.02302989, -0.020962052, -0.018608539,
        -0.015972615, -0.013059946, -0.009879752, -0.006443637, -0.002764011,
        0.00114632, 0.005272321, 0.009597827, 0.014107461, 0.018786028,
        0.023618005, 0.028587701, 0.033679874, 0.038879969, 0.044174387,
        0.04955013, 0.054994794, 0.060496516, 0.066043998, 0.071626168,
        0.077232297, 0.082852388, 0.088477845, 0.09410105, 0.099715411,
        0.105314674, 0.110892819, 0.116444212, 0.121963649, 0.127446089,
        0.132886752, 0.138281271, 0.143625983, 0.148918067, 0.154155521,
        0.159336978, 0.164461469, 0.169528598, 0.174538172, 0.179490209,
        0.184384927, 0.189222825, 0.194004593, 0.198731045, 0.203403085,
        0.208021636, 0.21258768, 0.217102204, 0.221566172, 0.22598057
      )    
      
      s_spline <- s_spline_lookup_table [age - 5]
      
      sigma <- exp (-3.171582 + 0.144358 * log (age) + s_spline)

    }
    
  }
  
  if (parameter == "tlc") {
    
    if (sex == 1) {
      
      s_spline_lookup_table <- c (
        0.088394747, 0.082392489, 0.076411414, 0.07045988, 0.06455393,
        0.0587139, 0.052960172, 0.047312491, 0.041782405, 0.036375791,
        0.031098408, 0.025955852, 0.020950808, 0.016083471, 0.011353956,
        0.006762462, 0.002311324, -0.001994877, -0.006151449, -0.010153679,
        -0.013995873, -0.017671069, -0.021172221, -0.024492279, -0.027624013,
        -0.030559899, -0.033292392, -0.035813943, -0.038117153, -0.040194918,
        -0.042040167, -0.043645829, -0.04500399, -0.046104636, -0.046937429,
        -0.047492034, -0.047758571, -0.047728609, -0.047394002, -0.0467466,
        -0.045779801, -0.044493167, -0.042887806, -0.040964825, -0.038725951,
        -0.036176071, -0.033321072, -0.030166843, -0.02671952, -0.022986855,
        -0.018977242, -0.014699081, -0.010160661, -0.005369375, -0.000332163,
        0.004944038, 0.010452182, 0.016184011, 0.022130513, 0.028282667,
        0.034631499, 0.041168795, 0.047886924, 0.054778268, 0.061835126,
        0.069048049, 0.076405924, 0.083897573, 0.091511779, 0.099236131,
        0.107056816, 0.114959945, 0.122931647, 0.130958981, 0.139030355
      )
      
      s_spline <- s_spline_lookup_table [age - 5]
      
      sigma <- exp (
        -2.0616143 +
        -0.0008534 * age +
        s_spline
      )
      
    } else {
      
      s_spline_lookup_table <- c (
        0.146232484, 0.137270933, 0.128326185, 0.119409995, 0.110542753,
        0.101745355, 0.093038695, 0.084443184, 0.075977267, 0.067658892,
        0.059506007, 0.051536187, 0.043760411, 0.036184003, 0.028812107,
        0.021649863, 0.014702784, 0.007977419, 0.001480491, -0.004781274,
        -0.010801036, -0.016570717, -0.022081474, -0.027324454, -0.032290795,
        -0.036970248, -0.041349864, -0.045416381, -0.049156537, -0.052556785,
        -0.05560155, -0.058274375, -0.060558797, -0.062438422, -0.063899779,
        -0.064933474, -0.065530413, -0.065681504, -0.065379434, -0.064625517,
        -0.063423665, -0.061777793, -0.059691902, -0.057172009, -0.054226143,
        -0.050862423, -0.047088968, -0.042914403, -0.038349042, -0.033403546,
        -0.028088575, -0.022415029, -0.016397049, -0.010051102, -0.003393707,
        0.003558609, 0.010787196, 0.018268524, 0.025978379, 0.033892544,
        0.041987161, 0.050241432, 0.058636127, 0.06715203, 0.075769948,
        0.084472213, 0.093243644, 0.102069289, 0.110934195, 0.119824077,
        0.128728506, 0.137638443, 0.146544845, 0.155438825, 0.164316207
      )      
      
      s_spline <- s_spline_lookup_table [age - 5]      
      
      sigma <- exp (
        -2.0999321 +
        0.0001564 * age +
        s_spline
      )      
      
    }
    
  }
  
  return (sigma)
  
}

get_specificity <- function (predicted, observed, cutoff) {
  
  fp <- 0 # False positives
  tn <- 0 # True negatives
  
  for (i in 1:length (observed)) {
    
    if (observed [i] == 0 & predicted [i] >= cutoff) fp <- fp + 1
    
    if (observed [i] == 0 & predicted [i] < cutoff) tn <- tn + 1
    
  }
  
  specificity <- tn / (fp + tn)
  
  return (specificity)
  
}

get_specificity_ers_ats <- function (data) {
  
  true_negatives <- data %>% 
    filter (fvc_z_score >= -1.645 & tlc_z_score >= -1.645) %>%
    nrow ()
    
  false_positives <- data %>% 
    filter (fvc_z_score < -1.645 & tlc_z_score > -1.645) %>%
    nrow () 

  numerator <- true_negatives
  
  denominator <- true_negatives + false_positives
  
  specificity <- numerator / denominator
  
  return (specificity)
  
}

get_tests <- function (data, model, label) {
  
  results <- rep (NA, 50)
  
  tests_bootstrap <- rep (NA, 10)
  
  for (i in 1:50) {
    
    lower <- i / 100
    
    upper <- 1 - lower
    
    predicted <- predict (model, data, type = 'prob')[,2]
    observed <- recode (pull (data, label), "yes" = 1, "no" = 0)
  
    actual <- tibble (
      predicted = predicted,
      observed = observed
    )
    
    tests <- actual %>% 
      filter (predicted > lower & predicted < upper) %>% 
      nrow () %>% 
      '/' (nrow (data))
      
    results [i] <- tests
    
  }
  
  return (results)
  
}

get_tests_ers_ats <- function (data) {
  
  tests <- data %>%
    mutate (predicted = case_when (
      fvc_z_score < -1.645 ~ 1,
      TRUE ~ 0)
    ) %>%
    filter (predicted == 1) %>%
    nrow () %>%
    '/' (nrow (data))
  
  return (tests)  
    
}

get_true_positive <- function (model, data, label, threshold) {
  
  predicted <- predict (model, data, type = 'prob')[,2]
  observed <- recode (pull (data, label), "yes" = 1, "no" = 0)
  
  actual <- tibble (
    data = data,
    predicted = predicted,
    observed = observed
  )
    
  true_positive <- actual %>% 
    filter (predicted >= threshold & observed == 1) %>% 
    nrow () %>% 
    '/' (nrow (data))
    
  return (true_positive)  
  
}

get_true_positive_ers_ats <- function (data) {
  
  data <- data %>% 
    mutate (predicted = case_when (
      fvc_z_score < -1.645 ~ 1,
      TRUE ~ 0)
    ) %>% 
    mutate (observed = case_when (
      tlc_z_score < -1.645 ~ 1,
      TRUE ~ 0)
    )
  
  true_positive <- data %>% 
    filter (predicted == 1 & observed == 1) %>% 
    nrow () %>% 
    '/' (nrow (data))
  
  return (true_positive)
  
}
