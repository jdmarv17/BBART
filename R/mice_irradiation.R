#' mice_irradiation
#'
#' This data comes from a Center for Medical Countermeasures Against Radiation Center funded
#' project to develop radiation exposure biomarkers and diagnostics for use in directing medical
#' treatment following a radiological incident where the radiation dose an individual is exposed to
#' would be unknown. C57Bl/6J mice, male and female, at 8 weeks of age, were exposed to
#' total body plus whole thorax irradiation using a 137Cs γ-ray source, or as a control received
#' no irradiation with identical handling (sham irradiation). Respiratory measures were recorded
#' over time using whole body plethysmography. The goal of the data collection were to develop a
#' means of deciphering radiation exposure by the resulting pulmonary response through identifying
#' respiratory measures that are associated with radiation exposure and building models that can
#' accurately predict radiation using respiratory measures.
#'
#'
#' @format A list with two elements:
#' \describe{
#'   \item{mice_covars}{A list of dataframes for subjects m1-m175 for 90 time points. x1 = "Respiratory frequency",
#'   x2 = "Tidal volume", x3 = "Minute ventilation", x4 = "Peak inspiratory flow", x5 = "Peak expiratory flow".}
#'   \item{y}{A vector containing the radiation status of each subject. 0 = "sham irradiation", 1 = "irradiation".}
#' }
#' @source Williams, Jacqueline P et al. “Animal models for medical countermeasures to radiation exposure”. In: Radiation research 173.4 (2010), pp. 557–578.
"mice_irradiation"
