###########################################
##                                       ##
##  MAPPING THE FIELD LITERATURE REVIEW  ##
##                                       ##
###########################################
##                                       ##
## Literature review data has undergone  ##
## initial processing; this completes    ##
## data processing, specifically         ##
## excluding low-count ERIC descriptors  ##
## and then identifies optimal number of ##
## clusters (subject to human review)    ##
## and visualizes clusters.              ##
##                                       ##
###########################################


###########################################
#         Load Required Libraries         #
###########################################


library(readr)
library(FactoMineR)
library(factoextra)
library(CAinterprTools)
library(ClustOfVar)
library(cluster)
library(NbClust)
library(vegan)
library(dendextend)
library(corrplot)
library(missMDA)
library(ggbiplot)
library(dplyr)
library(purrr)
library(safejoin)
library(tidyr)
library(stringr)

###########################################
#            Define Functions             #
###########################################

# This function creates a numeric binary 1-0 dataframe from the dataframe that is fed to it
binary.text.translate <- function(the.frame, start.col) {
  # Loop through each column
  for (i in start.col:length(the.frame)) {
    # Identify if it is a "no" incidence and if so, substitute it with a 0
    # Otherwise it must be a present incidence, so substitute it with a 1
    the.frame[[i]] <- if_else(str_detect(the.frame[[i]], "no_"), 0, 1)
    #the.frame[[i]] <- if_else(the.frame[[i]] <= 0, 0, 1)
    # Convert the column to numeric
    the.frame[[i]] <- as.numeric(the.frame[[i]])
  }
  # Send back the converted dataframe
  return(the.frame)
}

binary.numeric.translate <- function(the.frame, start.col) {
  # Loop through each column
  for (i in start.col:length(the.frame)) {
    # Identify if it is a "no" incidence and if so, substitute it with a 0
    # Otherwise it must be a present incidence, so substitute it with a 1
    #the.frame[[i]] <- if_else(substr(the.frame[[i]], 1, 3) == "no_", 0, 1)
    the.frame[[i]] <- if_else(the.frame[[i]] <= 0, 0, 1)
    # Convert the column to numeric
    the.frame[[i]] <- as.numeric(the.frame[[i]])
  }
  # Send back the converted dataframe
  return(the.frame)
}

###########################################
#                 Load Data               #
###########################################


# Load in data from identified traditions and add "traditions_" to column heading
# for identification purposes
traditions.frame <- read_csv("data/mtf-traditions.csv", col_names = TRUE, show_col_types = FALSE) %>%
  na.omit() %>% setNames(paste0('traditions_', names(.))) %>% dplyr::rename(manuscriptID = traditions_manuscriptID)

# Load in data from identified purposes and add "purposes_" to column heading
# for identification purposes
purposes.frame <- read_csv("data/mtf-purposes.csv", col_names = TRUE, show_col_types = FALSE) %>%
  na.omit() %>% setNames(paste0('purposes_', names(.))) %>% dplyr::rename("manuscriptID" = 1)

# Load in data from identified genres and add "genres_" to column heading
# for identification purposes
genres.frame <- read_csv("data/mtf-genres.csv", col_names = TRUE, show_col_types = FALSE) %>%
  na.omit() %>% setNames(paste0('genres_', names(.))) %>% dplyr::rename("manuscriptID" = 1)

# Load in data from identified data sources and add "sources_" to column heading
# for identification purposes
sources.frame <- read_csv("data/mtf-sources.csv", col_names = TRUE, show_col_types = FALSE) %>%
  na.omit() %>% setNames(paste0('sources_', names(.))) %>% dplyr::rename("manuscriptID" = 1)

# Load in date from identified ERIC descriptors and add "descriptors_" to column heading
# for identification purposes

# These are the ERIC descriptors from the first batch of articles, and then
# convert it to numeric (0 or 1)
keyword1.frame <- read.csv("data/mtf-eric_descriptors-1.csv", header = TRUE) %>%
  binary.text.translate(., 2)

# These are the ERIC descriptors from the second batch of articles, and then
# process them so that they are readable and able to be analyzed in a standardized manner
keyword2.frame <- read.csv("data/mtf-eric_descriptors-2.csv", header = FALSE)
keyword2.frame <- as.data.frame(lapply(keyword2.frame, function(keyword2.frame) gsub(" ", "_", keyword2.frame)))
keyword2.frame <- as.data.frame(lapply(keyword2.frame, function(keyword2.frame) gsub("[(]", "", keyword2.frame)))
keyword2.frame <- as.data.frame(lapply(keyword2.frame, function(keyword2.frame) gsub("[)]", "", keyword2.frame)))

# Convert to a long format for easier processing, and add the manuscriptID key
keyword2.frame <- pivot_longer(keyword2.frame, 2:27) %>%
  select(V1, value) %>%
  dplyr::rename("manuscriptID" = "V1", "name" = "value")
keyword2.frame$name <- tolower(keyword2.frame$name)
# Convert back to wide format coding presence as 1 and absence as 0
keyword2.frame <- keyword2.frame %>% filter(name != "") %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = name, values_from = value, id_cols = manuscriptID) %>%
  replace(is.na(.), 0)

# Pull the two ERIC descriptors lists together into one dataframe and process it
# and convert any NAs to 0 because a missing value means it is not there
descriptors.frame <- safe_full_join(keyword1.frame, keyword2.frame, by = "manuscriptID", conflict = coalesce) %>%
  replace(is.na(.), 0)

####                                             ####
# Identify descriptors to exclude to reduce "noise" #
####                                             ####

# The first list to exclude are those that deal with research methods since
# we human-coded these in a way that makes more sense to the analysis
methods.list <- c("action_research", "discourse_analysis", "factor_analysis",
                  "formative_evaluation", "measures_individuals",
                  "multiple_regression_analysis", "nonparametric_statistics",
                  "participatory_research", "reliability", "research_design",
                  "summative_evaluation", "validity", "case_studies",
                  "mixed_methods_research", "likert_scales", "interviews",
                  "surveys", "grounded_theory", "questionnaires",
                  "observation", "qualitative_research", "longitudinal_studies",
                  "correlation", "attitude_measures", "delphi_technique",
                  "data_collection", "data_analysis", "attitude_change",
                  "focus_groups", "semi_structured_interviews", "student_surveys",
                  "statistical_analysis", "teacher_surveys", "least_squares_statistics",
                  "regression_statistics", "units_of_study", "ethnography",
                  "researchers", "online_surveys", "comparative_analysis")

# Exclude descriptors that were search terms or just provide context that is not necessary
exclude.list <- c("college_students", "education_courses", "graduate_students",
                  "graduate_study", "higher_education", "influence_of_technology",
                  "majors_students", "masters_programs", "undergraduate_students",
                  "undergraduate_study", "technology_uses_in_education", "teacher_education",
                  "preservice_teacher_education", "educational_technology", "technology_education",
                  "preservice_teachers")

program.metadata.list <- c("program_descriptions", "program_development", "teacher_education_programs",
                           "summer_programs", "program_effectiveness", "program_implementation",
                           "program_design")

course.metadata.list <- c("course_content", "course_descriptions", "course_objectives",
                          "course_organization", "objectives", "workshops", "independent_study")

descriptors.frame <- descriptors.frame %>% select(-all_of(exclude.list)) %>%
  select(-all_of(methods.list)) %>%
  select(-all_of(program.metadata.list)) %>%
  select((-all_of(course.metadata.list)))

descriptors.frame <- descriptors.frame %>% select(order(colnames(descriptors.frame)))




attitudes.beliefs.list <- c("attitudes", "beliefs", "teacher_attitudes", "student_attitudes",
                            "administrator_attitudes", "student_teacher_attitudes")
descriptors.frame <- descriptors.frame %>% rowwise() %>% 
  dplyr::mutate(attitudes_beliefs = sum(c_across(any_of(attitudes.beliefs.list)))) %>%
  select(-all_of(attitudes.beliefs.list))

educational.change.list <- c("educational_innovation", "educational_trends", "experimental_programs",
                             "instructional_innovation", "intervention", "transformational_leadership")
descriptors.frame <- descriptors.frame %>% rowwise() %>% 
  dplyr::mutate(education_change = sum(c_across(any_of(educational.change.list)))) %>%
  select(-all_of(educational.change.list))

ela.list <- c("childrens_literature", "comprehension", "content_area_reading",
              "content_area_writing", "english_instruction", "language_arts", "literacy",
              "literacy_education", "poetry", "reader_response", "reading_assignments",
              "reading_instruction", "reading_strategies", "writing_across_the_curriculum",
              "writing_assignments", "writing_composition", "writing_instruction",
              "writing_processes", "writing_skills", "story_telling", "reading_skills",
              "english_teachers", "english_instruction", "english")
descriptors.frame <- descriptors.frame %>% rowwise() %>% 
  dplyr::mutate(ela = sum(c_across(any_of(ela.list)))) %>%
  select(-all_of(ela.list))

field.experiences.list <- c("college_school_cooperation", "field_experience_programs",
                            "internship_programs", "partnerships_in_education",
                            "professional_development_schools", "student_teaching", "practicums",
                            "service_learning", "mentors")
descriptors.frame <- descriptors.frame %>% rowwise() %>% 
  dplyr::mutate(field_experiences = sum(c_across(any_of(field.experiences.list)))) %>%
  select(-all_of(field.experiences.list))

instructional.design.list <- c("assignments", "curriculum_development", "curriculum_implementation",
                               "educational_objectives", "educational_planning",
                               "instructional_material_evaluation", "learning_activities",
                               "lesson_plans", "planning", "assignments", "curriculum_development",
                               "curriculum_implementation", "educational_objectives",
                               "worksheets", "textbooks", "instructional_development")
descriptors.frame <- descriptors.frame %>% rowwise() %>% 
  dplyr::mutate(instructional_design = sum(c_across(any_of(instructional.design.list)))) %>%
  select(-all_of(instructional.design.list))

learning.systems.list <- c("computer_assisted_instruction", "courseware", "integrated_learning_systems",
                           "web_based_instruction", "online_courses", "electronic_learning")
descriptors.frame <- descriptors.frame %>% rowwise() %>% 
  dplyr::mutate(learning_systems = sum(c_across(any_of(learning.systems.list)))) %>%
  select(-all_of(learning.systems.list))

metacognition.list <- c("cognitive_development", "decision_making", "prior_learning",
                        "problem_solving", "thinking_skills", "transfer_of_training", "value_judgment",
                        "self_evaluation_individuals", "self_evaluation_groups",
                        "self_concept", "modeling_psychology", "role_models")
descriptors.frame <- descriptors.frame %>% rowwise() %>% 
  dplyr::mutate(metacognition = sum(c_across(any_of(metacognition.list)))) %>%
  select(-all_of(metacognition.list))

simulations.games.list <- c("computer_simulation", "games", "simulated_environment")
descriptors.frame <- descriptors.frame %>% rowwise() %>% 
  dplyr::mutate(simulations_games = sum(c_across(any_of(simulations.games.list)))) %>%
  select(-all_of(simulations.games.list))

specific.technologies.list <-c ("audio_equipment", "computer_mediated_communication",
                                "computer_software", "electronic_equipment", "electronic_publishing",
                                "handheld_devices", "internet", "laptop_computers",
                                "open_source_technology", "telecommunications",
                                "video_technology", "visual_aids",
                                "web_2.0_technologies", "web_sites", "usability",
                                "use_studies")
descriptors.frame <- descriptors.frame %>% rowwise() %>% 
  dplyr::mutate(specific_technologies = sum(c_across(any_of(specific.technologies.list)))) %>%
  select(-all_of(specific.technologies.list))

sped.list <- c("special_education_teachers", "inclusion", "individualized_instruction",
               "disabilities")
descriptors.frame <- descriptors.frame %>% rowwise() %>% 
  dplyr::mutate(sped = sum(c_across(any_of(sped.list)))) %>%
  select(-all_of(sped.list))

stem.education.list <- c("astronomy", "elementary_school_science", "inquiry", "investigations",
                         "science_activities", "science_education", "secondary_school_science",
                         "science_teachers", "science_instruction", "scientific_methodology",
                         "physics", "misconceptions", "mathematics_teachers", "mathematics_instruction",
                         "mathematics_curriculum", "elementary_school_mathematics",
                         "stem_education")
descriptors.frame <- descriptors.frame %>% rowwise() %>% 
  dplyr::mutate(stem_education_sum = sum(c_across(any_of(stem.education.list)))) %>%
  select(-all_of(stem.education.list)) %>% dplyr::rename(stem_education = stem_education_sum)

structural.technology.issues.list <- c("access_to_computers", "barriers", "financial_support",
                                       "troubleshooting")
descriptors.frame <- descriptors.frame %>% rowwise() %>% 
  dplyr::mutate(structural_technology_issues = sum(c_across(any_of(structural.technology.issues.list)))) %>%
  select(-all_of(structural.technology.issues.list))

teacher.educators.list <- c("college_faculty", "faculty", "faculty_development")
descriptors.frame <- descriptors.frame %>% rowwise() %>% 
  dplyr::mutate(teacher_educators = sum(c_across(any_of(teacher.educators.list)))) %>%
  select(-all_of(teacher.educators.list))

teacher.effectiveness.list <- c("achievement_gains", "alignment_education", "competence",
                                "competency_based_teacher_education", "expertise",
                                "instructional_effectiveness", "performance_factors",
                                "student_teacher_evaluation", "teacher_competencies",
                                "teacher_competency_testing", "evaluation_methods", "success",
                                "scores", "pretests_posttests", "portfolios_background_materials")
descriptors.frame <- descriptors.frame %>% rowwise() %>% 
  dplyr::mutate(teacher_effectiveness = sum(c_across(any_of(teacher.effectiveness.list)))) %>%
  select(-all_of(teacher.effectiveness.list))

technology.literacy.list <- c("computer_literacy", "information_literacy", "information_skills",
                              "media_literacy", "multiple_literacies", "technological_literacy",
                              "technology_transfer")
descriptors.frame <- descriptors.frame %>% rowwise() %>% 
  dplyr::mutate(technology_literacy = sum(c_across(any_of(technology.literacy.list)))) %>%
  select(-all_of(technology.literacy.list))

tpack.list <- c("tpack", "knowledge", "pedagogical_content_knowledge", "knowledge_level")
descriptors.frame <- descriptors.frame %>% rowwise() %>% 
  dplyr::mutate(tpack_sum = sum(c_across(any_of(tpack.list)))) %>%
  select(-all_of(tpack.list)) %>% dplyr::rename(tpack = tpack_sum)

descriptors.frame <- descriptors.frame[colSums(descriptors.frame >= 2) > 0]
descriptors.frame <- as.data.frame(descriptors.frame)
descriptors.frame <- descriptors.frame %>% setNames(paste0('descriptors_', names(.))) %>%
  dplyr::rename("manuscriptID" = "descriptors_manuscriptID")


# Create the data frame for analysis by joining the individual data frames together
analysis.frame <- binary.text.translate(sources.frame, 2) %>%
  right_join(binary.text.translate(traditions.frame, 2), by = "manuscriptID") %>%
  right_join(binary.text.translate(genres.frame, 2), by = "manuscriptID") %>%
  right_join(binary.text.translate(purposes.frame, 2), by = "manuscriptID") %>%
  right_join(descriptors.frame, by = "manuscriptID") %>%
  as.data.frame() %>% na.omit()

# Set the row names as the Manuscript IDs (AUTHOR0000) and then remove that column
# because it is there solely to identify the rows
row.names(analysis.frame) <- analysis.frame$manuscriptID
analysis.frame$manuscriptID <- NULL


# ...and create a scaled Jaccard similarity index for analysis
# This Jaccard similarity index is necessary for hierarchical clustering analysis
analysis.jaccard <- scale(vegdist(analysis.frame, method = "jaccard"))

# Create a dataframe from the *variables* rather than the articles for PCA review
# This is accomplished by transposing the dataframe and scaling the creation of a
# Jaccard similarity index
variables.jaccard <- scale(vegdist(t(analysis.frame), method = "jaccard"))


###########################################
#       Run Initial PCA Analysis          #
###########################################

#pdf("initial_pca-plot.pdf", width = 36, height = 36)
#  analysis.pca <- PCA(analysis.frame)
#dev.off()

# condense into scales
analysis.scales <- analysis.frame
analysis.scales$manuscriptID <- rownames(analysis.scales)
rownames(analysis.scales) <- NULL
analysis.scales <- analysis.scales %>% relocate(manuscriptID)
social.scale.list <- c("traditions_deliberative", "genres_interpretive", "purposes_socialization")
analysis.scales <- analysis.scales %>% rowwise() %>% 
  dplyr::mutate(scale_social = sum(c_across(any_of(social.scale.list)))) %>%
  select(-all_of(social.scale.list))

scale.curriculum.list <- c("purposes_academic", "traditions_academic", "genres_other",
                      "genres_design", "sources_review")
analysis.scales <- analysis.scales %>% rowwise() %>% 
  dplyr::mutate(scale_curriculum = sum(c_across(any_of(scale.curriculum.list)))) %>%
  select(-all_of(scale.curriculum.list))

scale.justice.list <- c("purposes_justice", "purposes_other_purpose", 
                        "traditions_reconstruction", "traditions_other_tradition")
analysis.scales <- analysis.scales %>% rowwise() %>% 
  dplyr::mutate(scale_justice = sum(c_across(any_of(scale.justice.list)))) %>%
  select(-all_of(scale.justice.list))
analysis.scales <- analysis.scales %>% select(-scale_justice)

scale.pipeline.list <- c("traditions_technical", "genres_effects", "descriptors_tpack",
                         "sources_survey", "purposes_capital", "descriptors_teacher_effectiveness",
                         "descriptors_metacognition", "traditions_developmentalist")
analysis.scales <- analysis.scales %>% rowwise() %>% 
  dplyr::mutate(scale_pipeline = sum(c_across(any_of(scale.pipeline.list)))) %>%
  select(-all_of(scale.pipeline.list))

scale.methods.list <- c("sources_observation", "descriptors_teacher_educators",
                        "sources_artifact")
analysis.scales <- analysis.scales %>% rowwise() %>% 
  dplyr::mutate(scale_methods = sum(c_across(any_of(scale.methods.list)))) %>%
  select(-all_of(scale.methods.list))

scale.phenomenology.list <- c("genres_practitioner", "sources_interview", "purposes_subjectification")
analysis.scales <- analysis.scales %>% rowwise() %>% 
  dplyr::mutate(scale_phenomenology = sum(c_across(any_of(scale.phenomenology.list)))) %>%
  select(-all_of(scale.phenomenology.list))

scale.tech.in.field.list <- c("descriptors_field_experiences", "descriptors_technology_literacy")
analysis.scales <- analysis.scales %>% rowwise() %>% 
  dplyr::mutate(scale_tech_in_field = sum(c_across(any_of(scale.tech.in.field.list)))) %>%
  select(-all_of(scale.tech.in.field.list))

#scale.sol.list <- c("descriptors_instructional_design", "descriptors_stem_education",
#                    "descriptors_sped", "scale_curriculum")
#analysis.scales <- analysis.scales %>% rowwise() %>% 
#  dplyr::mutate(scale_science_of_teaching_learning = sum(c_across(any_of(scale.sol.list)))) %>%
#  select(-all_of(scale.sol.list))

scale.subjects.list <- c("descriptors_stem_education", "descriptors_sped", "scale_curriculum",
                    "descriptors_ela")
analysis.scales <- analysis.scales %>% rowwise() %>% 
  dplyr::mutate(scale_subjects = sum(c_across(any_of(scale.subjects.list)))) %>%
  select(-all_of(scale.subjects.list))

scale.social.tech.list <- c("scale_social", "scale_tech_in_field")
analysis.scales <- analysis.scales %>% rowwise() %>% 
  dplyr::mutate(scale_social_tech = sum(c_across(any_of(scale.social.tech.list)))) %>%
  select(-all_of(scale.social.tech.list))

analysis.scales <- as.data.frame(analysis.scales)
rownames(analysis.scales) <- analysis.scales$manuscriptID
analysis.scales$manuscriptID <- NULL
analysis.scales <- scale(analysis.scales)

library(hrbrthemes)

iu.colors <- c("#990000", "#FFAA00", "#056E41", "#006298", "#59264D")

# Compute hierarchical clustering on principal components
colnames(analysis.scales) <- c("instructional_design", "specific_technologies",
                               "attitudes_beliefs", "simulations_games",
                               "teacher_pipeline", "methods", "phenomonology",
                               "subject_areas", "socializing_technology_use")

analysis.pca <- PCA(analysis.scales, graph = FALSE)
res.hcpc <- HCPC(analysis.pca, graph = FALSE)

fviz_dend(res.hcpc,
          cex = 0.7,                     # Label size
          palette = iu.colors,               # Color palette see ?ggpubr::ggpar
          rect = FALSE, rect_fill = FALSE, # Add rectangle around groups
          ggtheme = theme_minimal(),
          #rect_border = iu.colors,           # Rectangle color
          labels_track_height = 0.8,      # Augment the room for labels
          main = "Dendrogram of Calculated Clusters") +
  scale_x_continuous(breaks = NULL)
ggsave("output/cluster_dend.pdf", width = 22, height = 6, units = "in", dpi = 300)

fviz_cluster(res.hcpc,
             repel = TRUE,            # Avoid label overlapping
             show.clust.cent = TRUE, # Show cluster centers
             palette = iu.colors,         # Color palette see ?ggpubr::ggpar
             ggtheme = theme_minimal(),
             main = "Cluster Map of Calculated Clusters")
ggsave("output/cluster_map.pdf", width = 11, height = 8.5, units = "in", dpi = 300, bg = "#FFFFFF")

fviz_eig(analysis.pca, addlabels = TRUE, ylim = c(0, 20), ggtheme = theme_minimal())

var <- get_pca_var(analysis.pca)
# Coordinates
var$coord
# Cos2: quality on the factore map
head(var$cos2)
# Contributions to the principal components
head(var$contrib)
# Coordinates of variables
head(var$coord, 4)
library(ggcorrplot)
library(ggsci)

  #corrplot(var$cos2, is.corr=FALSE, tl.col="black", method = "color", col = COL2('PuOr'))
ggcorrplot(t(var$cos2), method = "circle", colors = c("#990000", "#59264D", "#FFAA00")) +
  labs(x = "Dimensions", y = "Factors", title = "Dimension Correlations") +
  theme_minimal()
ggsave("output/corrplot.pdf", width = 8.5, height = 11, units = "in", dpi = 300)

# Total cos2 of variables on Dim.1 and Dim.2
fviz_cos2(analysis.pca, choice = "var", axes = 1)
fviz_cos2(analysis.pca, choice = "var", axes = 2)

# Color by cos2 values: quality on the factor map
fviz_pca_var(analysis.pca, col.var = "cos2",
             gradient.cols = c("#59264D", "#FFAA00", "#990000"), 
             repel = TRUE # Avoid text overlapping
)
ggsave("output/pca.pdf", width = 11, height = 11, units = "in", dpi = 300)



# Paragons!
res.hcpc$desc.ind$para
# Clusters
res.hcpc$data.clust
# Contributions to Clusters
res.hcpc$desc.var$quanti
# Descriptions of Dimensions
res.hcpc$desc.axes

cluster.details.1 <- res.hcpc$desc.var$quanti$`1`
# Cluster 1: Learning to Teach with Games and Simulations
cluster.details.2 <- res.hcpc$desc.var$quanti$`2`
# Cluster 2: Building a TPACK-Informed Teacher Workforce
cluster.details.3 <- res.hcpc$desc.var$quanti$`3`
# Cluster 3: Understanding What Teachers Do with Technology
cluster.details.4 <- res.hcpc$desc.var$quanti$`4`
# Cluster 4: Shaping Beliefs and Attitudes around Technology Use through Socialization
cluster.details.5 <- res.hcpc$desc.var$quanti$`5`
# Understanding Impressions Around Technology Use in the Field

head(var$contrib, 4)
corrplot(var$contrib, is.corr=FALSE)
# Contributions of variables to PC1
fviz_contrib(analysis.pca, choice = "var", axes = 1, top = 10)
# Contributions of variables to PC2
fviz_contrib(analysis.pca, choice = "var", axes = 2, top = 10)

