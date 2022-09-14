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
library(ggplot2)
library(FactoMineR)
library(factoextra)
library(vegan)
library(dendextend)
library(scales)
library(ggbiplot)
library(dplyr)
library(purrr)
library(safejoin)
library(tidyr)
library(stringr)
library(ggcorrplot)
library(ggalluvial)
library(ggfittext)
library(tidyclust)
library(flexclust)
library(ggpubr)
library(sna)
library(info.centrality)
library(igraph)
library(glue)
library(ggraph)
library(gt)

# Set file locations
project.path <- here::here()
data.path <- glue("{project.path}/data/")
output.path <- glue("{project.path}/output/")

###########################################
#            Define Functions             #
###########################################

# This function creates a numeric binary 1-0 dataframe from the dataframe that is fed to it
binary.text.translate <- function(the.frame, start.col) {
  # Loop through each column
  for (i in start.col:length(the.frame)) {
    # Identify if it is a "no" incidence and if so, substitute it with a 0
    # Otherwise it must be a present incidence, so substitute it with a 1
    the.frame[[i]] <-
      if_else(str_detect(the.frame[[i]], "no_"), 0, 1)
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

frame.factor.details <- function(the.frame, clusterNum) {
  the.frame$cluster <- clusterNum
  the.frame$factors <- rownames(the.frame)
  rownames(the.frame) <- NULL
  return(the.frame)
}

four_quadrant <- function(x, col_text = "#FFF7F8") {
  # Source: https://stackoverflow.com/questions/57923246/making-a-a-four-quadrant-proportional-area-chart-in-r
  nx <- length(x)
  sqx <- sqrt(x)
  df <- data.frame(
    x = c(sqx[1], -sqx[2], -sqx[3], sqx[4]) / 2,
    y = c(sqx[1], sqx[2], -sqx[3], -sqx[4]) / 2,
    size = sqx,
    label = x
  )
  mm <- max(df$size) * 1.1
  
  ggplot(data = df, aes(
    x = x,
    y = y,
    width = size,
    height = size,
    group = factor(size)
  )) +
    geom_tile(fill = iu.4.colors) +
    geom_text(aes(label = label), col = col_text, size = 5) +
    geom_hline(aes(yintercept = 0), size = 0.8) +
    geom_vline(aes(xintercept = 0), size = 0.8) +
    coord_fixed() +
    theme_void() +
    theme(legend.position = "none")
}

ego.effective.size <- function(g, ego, ...) {
  egonet <- induced.subgraph(g, neighbors(g, ego, ...))
  n <- vcount(egonet)
  t <- ecount(egonet)
  return(n - (2 * t) / n)
}

isolate.cluster <- function(the.cluster, the.frame) {
  the.frame <- the.frame %>%
    left_join(analysis.cluster.frame, by = "manuscriptID") %>%
    filter(cluster == the.cluster) %>%
    select(manuscriptID, cite)
  return(the.frame)
}

calculate.centrality <- function(the.frame) {
  
  the.graph <- simplify(
    graph_from_edgelist(as.matrix(the.frame), directed = FALSE),
    remove.loops = TRUE,
    remove.multiple = TRUE
  )
  
  data.degree <-
    centralization.degree(the.graph, mode = "all", loops = FALSE)
  data.net.degree <- the.graph$centralization
  data.betweenness <-
    centralization.betweenness(the.graph, directed = FALSE)
  data.net.betweenness <- data.betweenness$centralization
  data.closeness <-
    centralization.closeness(the.graph, mode = "all")
  data.net.closeness <- data.closeness$centralization
  data.eigen <- centralization.evcent(the.graph, directed = FALSE)
  data.coreness <- coreness(the.graph, mode = "all")
  data.info <- info.centrality.vertex(the.graph)
  
  analysis.network.data <- data.frame(
    degree = data.degree$res,
    betweenness = data.betweenness$res,
    closeness = data.closeness$res,
    eigen = data.eigen$vector,
    coreness = data.coreness,
    info = data.info
  )
  return(analysis.network.data)
}

calculate.keyactors <- function(the.frame, the.cluster) {
  the.filename <- glue("output/keyactors{the.cluster}-plot.pdf")
  key.frame <- the.frame %>%
    select(betweenness, eigen)
  key.frame$scaled_betweenness <-
    rescale(key.frame$betweenness, to = c(0.01, 0.99))
  key.frame$scaled_eigen <- rescale(key.frame$eigen, to = c(0, 1))
  key.frame$interaction_score <-
    (key.frame$scaled_betweenness * key.frame$scaled_eigen)
  key.frame$interaction_score <-
    rescale(key.frame$interaction_score, to = c(0, 1))
  key.frame$manuscriptID <- rownames(key.frame)
  rownames(key.frame) <- NULL
  key.frame$cluster <- the.cluster
  key.res <- lm(eigen ~ betweenness, data = key.frame)$residuals
  key.frame <- transform(key.frame, residuals = key.res)
  key.frame <- key.frame %>%
    select(
      manuscriptID,
      cluster,
      betweenness,
      scaled_betweenness,
      eigen,
      scaled_eigen,
      interaction_score,
      residuals
    ) %>%
    filter(!manuscriptID %in% analysis.cluster.frame$manuscriptID)
  key.ymedian <- median(key.frame$eigen)
  key.xmedian <- median(key.frame$scaled_betweenness)
  key.frame <- key.frame %>%
    mutate(keystatus = case_when((scaled_eigen > key.ymedian &
                                    scaled_betweenness > key.xmedian) ~ "Hub",
                                 (scaled_eigen > key.ymedian &
                                    scaled_betweenness < key.xmedian) ~ "Pulse-Taker",
                                 (scaled_eigen < key.ymedian &
                                    scaled_betweenness > key.ymedian) ~ "Gatekeeper"
    )) %>%
    na.omit()
  
  key.xmin <- min(key.frame$scaled_betweenness) - 0.1
  key.xmax <- max(key.frame$scaled_betweenness) + 0.1
  key.ymin <- min(key.frame$scaled_eigen) - 0.1
  key.ymax <- max(key.frame$scaled_eigen) + 0.1
  
  key.plot <-
    ggscatter(
      key.frame,
      x = "scaled_betweenness",
      y = "scaled_eigen",
      label = "manuscriptID",
      label.rectangle = FALSE,
      repel = TRUE,
      theme = theme_minimal(),
      xlab = "Betweenness Centrality (scaled)",
      ylab = "Eigenvector Centrality (scaled)",
      point = TRUE,
      ylim = c(key.ymin, key.ymax),
      xlim = c(key.xmin, key.xmax),
      show.legend.text = FALSE,
      color = "residuals",
      legend = "bottom",
      conf.int = FALSE, 
      cor.coef = FALSE
    ) +
    geom_smooth(method = "lm", color = "#EEEEEE", se = FALSE, fullrange = TRUE, linetype = "solid") +
    scale_color_gradient2(low = iu.gradient$low[the.cluster],
                          high = iu.gradient$high[the.cluster],
                          mid = iu.gradient$mid[the.cluster]) +
    geom_hline(yintercept = key.ymedian) +
    geom_vline(xintercept = key.xmedian) +
    geom_label(
      aes(
        x = key.xmax,
        y = key.ymin,
        label = "Gatekeepers",
        hjust = 1
      ),
      color = "#330D2B",
      fill = "#DECADC"
    ) +
    geom_label(
      aes(
        x = key.xmax,
        y = key.ymax,
        label = "Hubs",
        hjust = 1
      ),
      color = "#330D2B",
      fill = "#DECADC"
    )
  ggsave(
    key.plot,
    filename = the.filename,
    width = 11.5,
    height = 8,
    units = "in",
    dpi = 300
  )
  key.frame <- key.frame %>%
    #dplyr::group_by(keystatus) %>%
    arrange(desc(keystatus), desc(interaction_score), manuscriptID)# %>%
  #slice(1:20) %>% ungroup()
  return(key.frame)
}

network.plot <-
  function(the.frame,
           the.ka,
           the.cluster,
           the.size,
           the.slice,
           the.width,
           the.height) {
    the.filename <- glue("output/cites-network-{the.cluster}.pdf")
    #the.frame <- cites.2.frame
    #the.ka <- cites.2.keyactors
    #the.slice <- 20
    #the.cluster <- 1
    #the.size <- 4
    the.ka$interaction_score <-
      (rescale(the.ka$interaction_score, to = c(1, 15))) + the.size
    the.ka <- the.ka %>%
      #    arrange(desc(keystatus), desc(interaction_score), manuscriptID) %>%
      dplyr::group_by(keystatus) %>%
      slice(1:the.slice) %>% ungroup() %>%
      mutate(
        point_shape = case_when(
          keystatus == "Hub" ~ "circle",
          keystatus == "Gatekeeper" ~ "diamond",
          keystatus == "Pulse-Taker" ~ "triangle"
        )
      ) %>%
      filter(keystatus != "Pulse-Taker")
    
    the.frame <- the.frame %>%
      filter(cite %in% the.ka$manuscriptID)
    the.graph <- simplify(
      graph_from_edgelist(as.matrix(the.frame[, 1:2]), directed = FALSE),
      remove.loops = TRUE,
      remove.multiple = TRUE
    )
    
    V(the.graph)$color <-
      as.character(analysis.cluster.frame$color[match(V(the.graph)$name, analysis.cluster.frame$manuscriptID)])
    V(the.graph)$color[is.na(V(the.graph)$color)] <- "#A7A9AB"
    V(the.graph)$shape <-
      as.character(the.ka$point_shape[match(V(the.graph)$name, the.ka$manuscriptID)])
    V(the.graph)$shape[is.na(V(the.graph)$shape)] <- "square"
    V(the.graph)$interaction_score <-
      as.character(the.ka$interaction_score[match(V(the.graph)$name, the.ka$manuscriptID)])
    V(the.graph)$interaction_score[is.na(V(the.graph)$interaction_score)] <-
      as.character(mean(the.ka$interaction_score))
    #V(the.graph)$interaction_score <- rescale(as.numeric(V(the.graph)$interaction_score), to = c(4, 20))
    set.seed(123)
    the.plot <- ggraph(the.graph, layout = "graphopt") +
      geom_edge_link(color = "#EEEEEE") +
      geom_node_point(
        color = V(the.graph)$color,
        alpha = 0.9,
        size = as.numeric(V(the.graph)$interaction_score),
        shape = V(the.graph)$shape
      ) +
      geom_node_text(aes(label = name),
                     size = 3,
                     color = "#243142",
                     repel = TRUE) +
      theme_void()
    ggsave(
      the.filename,
      plot = the.plot,
      width = the.width,
      height = the.height,
      units = "in",
      dpi = 300
    )
  }

#########################
# Define Color Palettes #
#########################

# Define colors, using IU branding colors https://www.iu.edu/brand/brand-expression/visual-language/color/index.html
iu.colors <-
  c("#990000", "#FFAA00", "#056E41", "#006298", "#59264D")
iu.colors.fade <<-
  c("#990000", "#FFAA00", "#056E41", "#006298", "#59264D")
#iu.colors.fade <<- c("#F23A3F", "#FFD563", "#63B363", "#63B1D3", "#A07498")
#iu.4.colors <- c("#0062986C", "#59264D6C", "#9900006C", "#056E416C")
iu.4.colors <- c("#006298", "#59264D", "#990000", "#056E41")
iu.10.colors <-
  c(
    "#990000",
    "#FFAA00",
    "#056E41",
    "#A7A9AB",
    "#006298",
    "#59264D",
    "#DF3603",
    "#01426A",
    "#49AFC7"
  )
iu.crimson <-
  c(
    "#5A0C0C",
    "#800000",
    "#990000",
    "#DC231E",
    "#F23A3F",
    "#FF636A",
    "#FFD6DB",
    "#FFF7F8"
  )
iu.crimson <-
  c(
    "#FFF7F8",
    "#FFD6DB",
    "#FF636A",
    "#F23A3F",
    "#DC231E",
    "#990000",
    "#800000",
    "#5A0C0C"
  )
iu.gradient <<-
  data.frame(
    high = c(
      "#5A0C0C",
      "#613D00",
      "#004421",
      "#00385F",
      "#330D2B",
      "#4A3C31"
    ),
    mid = c(
      "#990000",
      "#FFAA00",
      "#056E41",
      "#006298",
      "#59264D",
      "#A7A9AB"
    ),
    low = c(
      "#FFD6DB",
      "#FFF4C6",
      "#DEE8C6",
      "#C6ECF6",
      "#DECADC",
      "#EEEEEE"
    )
  )


###########################################
#                 Load Data               #
###########################################


# Load in data from identified traditions and add "traditions_" to column heading
# for identification purposes
traditions.frame <-
  read_csv("data/mtf-traditions.csv",
           col_names = TRUE,
           show_col_types = FALSE) %>%
  na.omit() %>% setNames(paste0('traditions_', names(.))) %>% dplyr::rename(manuscriptID = traditions_manuscriptID)

# Load in data from identified purposes and add "purposes_" to column heading
# for identification purposes
purposes.frame <-
  read_csv("data/mtf-purposes.csv",
           col_names = TRUE,
           show_col_types = FALSE) %>%
  na.omit() %>% setNames(paste0('purposes_', names(.))) %>% dplyr::rename("manuscriptID" = 1)

# Load in data from identified genres and add "genres_" to column heading
# for identification purposes
genres.frame <-
  read_csv("data/mtf-genres.csv",
           col_names = TRUE,
           show_col_types = FALSE) %>%
  na.omit() %>% setNames(paste0('genres_', names(.))) %>% dplyr::rename("manuscriptID" = 1)

# Load in data from identified data sources and add "sources_" to column heading
# for identification purposes
sources.frame <-
  read_csv("data/mtf-sources.csv",
           col_names = TRUE,
           show_col_types = FALSE) %>%
  na.omit() %>% setNames(paste0('sources_', names(.))) %>% dplyr::rename("manuscriptID" = 1)

# Load in date from identified ERIC descriptors and add "descriptors_" to column heading
# for identification purposes

# These are the ERIC descriptors from the first batch of articles, and then
# convert it to numeric (0 or 1)
keyword1.frame <-
  read.csv("data/mtf-eric_descriptors-1.csv", header = TRUE) %>%
  binary.text.translate(., 2)

# These are the ERIC descriptors from the second batch of articles, and then
# process them so that they are readable and able to be analyzed in a standardized manner
keyword2.frame <-
  read.csv("data/mtf-eric_descriptors-2.csv", header = FALSE)
keyword2.frame <-
  as.data.frame(lapply(keyword2.frame, function(keyword2.frame)
    gsub(" ", "_", keyword2.frame)))
keyword2.frame <-
  as.data.frame(lapply(keyword2.frame, function(keyword2.frame)
    gsub("[(]", "", keyword2.frame)))
keyword2.frame <-
  as.data.frame(lapply(keyword2.frame, function(keyword2.frame)
    gsub("[)]", "", keyword2.frame)))

# Convert to a long format for easier processing, and add the manuscriptID key
keyword2.frame <- pivot_longer(keyword2.frame, 2:27) %>%
  select(V1, value) %>%
  dplyr::rename("manuscriptID" = "V1", "name" = "value")
keyword2.frame$name <- tolower(keyword2.frame$name)
# Convert back to wide format coding presence as 1 and absence as 0
keyword2.frame <- keyword2.frame %>% filter(name != "") %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = name,
              values_from = value,
              id_cols = manuscriptID) %>%
  replace(is.na(.), 0)

# Pull the two ERIC descriptors lists together into one dataframe and process it
# and convert any NAs to 0 because a missing value means it is not there
descriptors.frame <-
  safe_full_join(keyword1.frame,
                 keyword2.frame,
                 by = "manuscriptID",
                 conflict = coalesce) %>%
  replace(is.na(.), 0)

#####################################################
# Identify descriptors to exclude to reduce "noise" #
#####################################################

# The first list to exclude are those that deal with research methods since
# we human-coded these in a way that makes more sense to the analysis
methods.list <-
  c(
    "action_research",
    "discourse_analysis",
    "factor_analysis",
    "formative_evaluation",
    "measures_individuals",
    "multiple_regression_analysis",
    "nonparametric_statistics",
    "participatory_research",
    "reliability",
    "research_design",
    "summative_evaluation",
    "validity",
    "case_studies",
    "mixed_methods_research",
    "likert_scales",
    "interviews",
    "surveys",
    "grounded_theory",
    "questionnaires",
    "observation",
    "qualitative_research",
    "longitudinal_studies",
    "correlation",
    "attitude_measures",
    "delphi_technique",
    "data_collection",
    "data_analysis",
    "attitude_change",
    "focus_groups",
    "semi_structured_interviews",
    "student_surveys",
    "statistical_analysis",
    "teacher_surveys",
    "least_squares_statistics",
    "regression_statistics",
    "units_of_study",
    "ethnography",
    "researchers",
    "online_surveys",
    "comparative_analysis"
  )

# Exclude descriptors that were search terms or just provide context that is not necessary
exclude.list <-
  c(
    "college_students",
    "education_courses",
    "graduate_students",
    "graduate_study",
    "higher_education",
    "influence_of_technology",
    "majors_students",
    "masters_programs",
    "undergraduate_students",
    "undergraduate_study",
    "technology_uses_in_education",
    "teacher_education",
    "preservice_teacher_education",
    "educational_technology",
    "technology_education",
    "preservice_teachers"
  )

# Exclude extraneous descriptors that address programmatic and course items
program.metadata.list <-
  c(
    "program_descriptions",
    "program_development",
    "teacher_education_programs",
    "summer_programs",
    "program_effectiveness",
    "program_implementation",
    "program_design"
  )

course.metadata.list <-
  c(
    "course_content",
    "course_descriptions",
    "course_objectives",
    "course_organization",
    "objectives",
    "workshops",
    "independent_study"
  )

# Actually remove the descriptors here...
descriptors.frame <-
  descriptors.frame %>% select(-all_of(exclude.list)) %>%
  select(-all_of(methods.list)) %>%
  select(-all_of(program.metadata.list)) %>%
  select((-all_of(course.metadata.list)))

# ...And order them alphabetically
descriptors.frame <-
  descriptors.frame %>% select(order(colnames(descriptors.frame)))

######################################################
#         Initial Combination of Descriptors         #
######################################################

# This reduction of descriptors involves combining original descriptors into categories.
# The categories are summed so that a sense of scale is retained. This also allows
# them to be continuous variables for PCA analysis. Even though some descriptors can
# be matched to specific fields and efforts in education, the definitions are considered
# broad rather than narrow.

# Combine descriptors having to do with attitudes and beliefs
attitudes.beliefs.list <-
  c(
    "attitudes",
    "beliefs",
    "teacher_attitudes",
    "student_attitudes",
    "administrator_attitudes",
    "student_teacher_attitudes"
  )
descriptors.frame <- descriptors.frame %>% rowwise() %>%
  dplyr::mutate(attitudes_beliefs = sum(c_across(any_of(
    attitudes.beliefs.list
  )))) %>%
  select(-all_of(attitudes.beliefs.list))

# Combine descriptors having to do with educational change efforts
educational.change.list <-
  c(
    "educational_innovation",
    "educational_trends",
    "experimental_programs",
    "instructional_innovation",
    "intervention",
    "transformational_leadership"
  )
descriptors.frame <- descriptors.frame %>% rowwise() %>%
  dplyr::mutate(education_change = sum(c_across(any_of(
    educational.change.list
  )))) %>%
  select(-all_of(educational.change.list))

# Combine descriptors having to do with English-Language Arts
ela.list <-
  c(
    "childrens_literature",
    "comprehension",
    "content_area_reading",
    "content_area_writing",
    "english_instruction",
    "language_arts",
    "literacy",
    "literacy_education",
    "poetry",
    "reader_response",
    "reading_assignments",
    "reading_instruction",
    "reading_strategies",
    "writing_across_the_curriculum",
    "writing_assignments",
    "writing_composition",
    "writing_instruction",
    "writing_processes",
    "writing_skills",
    "story_telling",
    "reading_skills",
    "english_teachers",
    "english_instruction",
    "english"
  )
descriptors.frame <- descriptors.frame %>% rowwise() %>%
  dplyr::mutate(ela = sum(c_across(any_of(ela.list)))) %>%
  select(-all_of(ela.list))

# Combine descriptors having to do with field experiences specifically
field.experiences.list <-
  c(
    "college_school_cooperation",
    "field_experience_programs",
    "internship_programs",
    "partnerships_in_education",
    "professional_development_schools",
    "student_teaching",
    "practicums",
    "service_learning",
    "mentors"
  )
descriptors.frame <- descriptors.frame %>% rowwise() %>%
  dplyr::mutate(field_experiences = sum(c_across(any_of(
    field.experiences.list
  )))) %>%
  select(-all_of(field.experiences.list))

# Combine descriptors having to do with Instructional Design
instructional.design.list <-
  c(
    "assignments",
    "curriculum_development",
    "curriculum_implementation",
    "educational_objectives",
    "educational_planning",
    "instructional_material_evaluation",
    "learning_activities",
    "lesson_plans",
    "planning",
    "assignments",
    "curriculum_development",
    "curriculum_implementation",
    "educational_objectives",
    "worksheets",
    "textbooks",
    "instructional_development"
  )
descriptors.frame <- descriptors.frame %>% rowwise() %>%
  dplyr::mutate(instructional_design = sum(c_across(any_of(
    instructional.design.list
  )))) %>%
  select(-all_of(instructional.design.list))

# Combine descriptors having to do with learing systems, online learning, etc.
learning.systems.list <-
  c(
    "computer_assisted_instruction",
    "courseware",
    "integrated_learning_systems",
    "web_based_instruction",
    "online_courses",
    "electronic_learning"
  )
descriptors.frame <- descriptors.frame %>% rowwise() %>%
  dplyr::mutate(learning_systems = sum(c_across(any_of(
    learning.systems.list
  )))) %>%
  select(-all_of(learning.systems.list))

# Combine descriptors having to do with metacognition and psychology
metacognition.list <-
  c(
    "cognitive_development",
    "decision_making",
    "prior_learning",
    "problem_solving",
    "thinking_skills",
    "transfer_of_training",
    "value_judgment",
    "self_evaluation_individuals",
    "self_evaluation_groups",
    "self_concept",
    "modeling_psychology",
    "role_models"
  )
descriptors.frame <- descriptors.frame %>% rowwise() %>%
  dplyr::mutate(metacognition = sum(c_across(any_of(metacognition.list)))) %>%
  select(-all_of(metacognition.list))

# Combine descriptors having to do with simulations and games
simulations.games.list <-
  c("computer_simulation", "games", "simulated_environment")
descriptors.frame <- descriptors.frame %>% rowwise() %>%
  dplyr::mutate(simulations_games = sum(c_across(any_of(
    simulations.games.list
  )))) %>%
  select(-all_of(simulations.games.list))

# Combine descriptors having to do with specific technologies and how they are being used/evaluated
specific.technologies.list <-
  c (
    "audio_equipment",
    "computer_mediated_communication",
    "computer_software",
    "electronic_equipment",
    "electronic_publishing",
    "handheld_devices",
    "internet",
    "laptop_computers",
    "open_source_technology",
    "telecommunications",
    "video_technology",
    "visual_aids",
    "web_2.0_technologies",
    "web_sites",
    "usability",
    "use_studies"
  )
descriptors.frame <- descriptors.frame %>% rowwise() %>%
  dplyr::mutate(specific_technologies = sum(c_across(any_of(
    specific.technologies.list
  )))) %>%
  select(-all_of(specific.technologies.list))

# Combine descriptors having to do with special education
sped.list <-
  c(
    "special_education_teachers",
    "inclusion",
    "individualized_instruction",
    "disabilities"
  )
descriptors.frame <- descriptors.frame %>% rowwise() %>%
  dplyr::mutate(sped = sum(c_across(any_of(sped.list)))) %>%
  select(-all_of(sped.list))

# Combine descriptors having to do with science, technology, engineering, and mathematics education
stem.education.list <-
  c(
    "astronomy",
    "elementary_school_science",
    "inquiry",
    "investigations",
    "science_activities",
    "science_education",
    "secondary_school_science",
    "science_teachers",
    "science_instruction",
    "scientific_methodology",
    "physics",
    "misconceptions",
    "mathematics_teachers",
    "mathematics_instruction",
    "mathematics_curriculum",
    "elementary_school_mathematics",
    "stem_education"
  )
descriptors.frame <- descriptors.frame %>% rowwise() %>%
  dplyr::mutate(stem_education_sum = sum(c_across(any_of(
    stem.education.list
  )))) %>%
  select(-all_of(stem.education.list)) %>% dplyr::rename(stem_education = stem_education_sum)

# Combine descriptors having to do with issues that are structural in nature
structural.technology.issues.list <-
  c("access_to_computers",
    "barriers",
    "financial_support",
    "troubleshooting")
descriptors.frame <- descriptors.frame %>% rowwise() %>%
  dplyr::mutate(structural_technology_issues = sum(c_across(
    any_of(structural.technology.issues.list)
  ))) %>%
  select(-all_of(structural.technology.issues.list))

# Combine descriptors having to do with teacher educators, not just teacher education
teacher.educators.list <-
  c("college_faculty", "faculty", "faculty_development")
descriptors.frame <- descriptors.frame %>% rowwise() %>%
  dplyr::mutate(teacher_educators = sum(c_across(any_of(
    teacher.educators.list
  )))) %>%
  select(-all_of(teacher.educators.list))

# Combine descriptors having to do with teacher effectiveness
teacher.effectiveness.list <-
  c(
    "achievement_gains",
    "alignment_education",
    "competence",
    "competency_based_teacher_education",
    "expertise",
    "instructional_effectiveness",
    "performance_factors",
    "student_teacher_evaluation",
    "teacher_competencies",
    "teacher_competency_testing",
    "evaluation_methods",
    "success",
    "scores",
    "pretests_posttests",
    "portfolios_background_materials"
  )
descriptors.frame <- descriptors.frame %>% rowwise() %>%
  dplyr::mutate(teacher_effectiveness = sum(c_across(any_of(
    teacher.effectiveness.list
  )))) %>%
  select(-all_of(teacher.effectiveness.list))

# Combine descriptors having to do with technology literacy
technology.literacy.list <-
  c(
    "computer_literacy",
    "information_literacy",
    "information_skills",
    "media_literacy",
    "multiple_literacies",
    "technological_literacy",
    "technology_transfer"
  )
descriptors.frame <- descriptors.frame %>% rowwise() %>%
  dplyr::mutate(technology_literacy = sum(c_across(any_of(
    technology.literacy.list
  )))) %>%
  select(-all_of(technology.literacy.list))

# Combine descriptors having to do with TPACK, a well-known framework in the field
tpack.list <-
  c("tpack",
    "knowledge",
    "pedagogical_content_knowledge",
    "knowledge_level")
descriptors.frame <- descriptors.frame %>% rowwise() %>%
  dplyr::mutate(tpack_sum = sum(c_across(any_of(tpack.list)))) %>%
  select(-all_of(tpack.list)) %>% dplyr::rename(tpack = tpack_sum)

# Reduce descriptors by removing any that apply to two or fewer articles.
# There is a fine line between enough data and too much data, and this step is
# done to set the balance to enough data without losing a sense of how the
# articles connect with one another.
descriptors.frame <-
  descriptors.frame[colSums(descriptors.frame >= 2) > 0]
descriptors.frame <- as.data.frame(descriptors.frame)
# Add "descriptors_" tag to column names to identify them as descriptors
descriptors.frame <-
  descriptors.frame %>% setNames(paste0('descriptors_', names(.))) %>%
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
analysis.jaccard <-
  scale(vegdist(analysis.frame, method = "jaccard"))

# Create a dataframe from the *variables* rather than the articles for PCA review
# This is accomplished by transposing the dataframe and scaling the creation of a
# Jaccard similarity index
variables.jaccard <-
  scale(vegdist(t(analysis.frame), method = "jaccard"))


###########################################
#       Run Initial PCA Analysis          #
###########################################

# An initial Primary Component Analysis is conducted to identify overlapping
# factors for further reduction based on the visualization.

initial.analysis.pca <- PCA(analysis.frame, graph = FALSE)
fviz_pca_var(
  initial.analysis.pca,
  col.var = "cos2",
  gradient.cols = c("#59264D", "#FFAA00", "#990000"),
  repel = TRUE,
  # Avoid text overlapping
)
ggsave(
  "output/initial-pca.pdf",
  width = 8.5,
  height = 8.5,
  units = "in",
  dpi = 300
)

###########################################
#    Further Condense Data Into Scales    #
###########################################
# Based on examination of the PCA visualization, factors from across the full
# dataframe are reduced by combining them into scales. Like the process with the
# descriptors only, these are summed to provide a sense of scale and therefore
# represent continuous variables rather than categorical ones.

# Prepare the dataframe for scaling by first making a copy
analysis.scales <- analysis.frame
# A manuscriptID key column is created because column names are removed during
# dplyr mutate procedures and then moved to the first position.
analysis.scales$manuscriptID <- rownames(analysis.scales)
rownames(analysis.scales) <- NULL
analysis.scales <- analysis.scales %>% relocate(manuscriptID)

# Combine factors that connect with socialization into the teaching field
social.scale.list <-
  c("traditions_deliberative",
    "genres_interpretive",
    "purposes_socialization")
analysis.scales <- analysis.scales %>% rowwise() %>%
  dplyr::mutate(scale_social = sum(c_across(any_of(social.scale.list)))) %>%
  select(-all_of(social.scale.list))

# Combine factors that connect with curriculum and academics
scale.curriculum.list <-
  c(
    "purposes_academic",
    "traditions_academic",
    "genres_other",
    "genres_design",
    "sources_review"
  )
analysis.scales <- analysis.scales %>% rowwise() %>%
  dplyr::mutate(scale_curriculum = sum(c_across(any_of(
    scale.curriculum.list
  )))) %>%
  select(-all_of(scale.curriculum.list))

# Combine factors that connect with social justice
scale.justice.list <-
  c(
    "purposes_justice",
    "purposes_other_purpose",
    "traditions_reconstruction",
    "traditions_other_tradition"
  )
analysis.scales <- analysis.scales %>% rowwise() %>%
  dplyr::mutate(scale_justice = sum(c_across(any_of(scale.justice.list)))) %>%
  select(-all_of(scale.justice.list))
analysis.scales <- analysis.scales %>% select(-scale_justice)

# Combine factors that connect with the "Teacher Preparation Pipeline"
scale.pipeline.list <-
  c(
    "traditions_technical",
    "genres_effects",
    "descriptors_tpack",
    "sources_survey",
    "purposes_capital",
    "descriptors_teacher_effectiveness",
    "descriptors_metacognition",
    "traditions_developmentalist"
  )
analysis.scales <- analysis.scales %>% rowwise() %>%
  dplyr::mutate(scale_pipeline = sum(c_across(any_of(
    scale.pipeline.list
  )))) %>%
  select(-all_of(scale.pipeline.list))

# Combine factors that connect with methods courses
scale.methods.list <-
  c("sources_observation",
    "descriptors_teacher_educators",
    "sources_artifact")
analysis.scales <- analysis.scales %>% rowwise() %>%
  dplyr::mutate(scale_methods = sum(c_across(any_of(scale.methods.list)))) %>%
  select(-all_of(scale.methods.list))

# Combine factors that connect with understanding with a phenomological perspective
scale.phenomenology.list <-
  c("genres_practitioner",
    "sources_interview",
    "purposes_subjectification")
analysis.scales <- analysis.scales %>% rowwise() %>%
  dplyr::mutate(scale_phenomenology = sum(c_across(any_of(
    scale.phenomenology.list
  )))) %>%
  select(-all_of(scale.phenomenology.list))

# Combine factors that connect with knowing how to use technology in authentic experiences
scale.tech.in.field.list <-
  c("descriptors_field_experiences",
    "descriptors_technology_literacy")
analysis.scales <- analysis.scales %>% rowwise() %>%
  dplyr::mutate(scale_tech_in_field = sum(c_across(any_of(
    scale.tech.in.field.list
  )))) %>%
  select(-all_of(scale.tech.in.field.list))

# Combine factors that connect with specific academic subjects
scale.subjects.list <-
  c(
    "descriptors_stem_education",
    "descriptors_sped",
    "scale_curriculum",
    "descriptors_ela"
  )
analysis.scales <- analysis.scales %>% rowwise() %>%
  dplyr::mutate(scale_subjects = sum(c_across(any_of(
    scale.subjects.list
  )))) %>%
  select(-all_of(scale.subjects.list))

# Combine factors that connect with socializing into authentic technology use
scale.social.tech.list <- c("scale_social", "scale_tech_in_field")
analysis.scales <- analysis.scales %>% rowwise() %>%
  dplyr::mutate(scale_social_tech = sum(c_across(any_of(
    scale.social.tech.list
  )))) %>%
  select(-all_of(scale.social.tech.list))

list.of.lists <-
  do.call(
    c,
    list(
      tpack.list,
      technology.literacy.list,
      teacher.effectiveness.list,
      teacher.educators.list,
      structural.technology.issues.list,
      stem.education.list,
      sped.list,
      specific.technologies.list,
      simulations.games.list,
      metacognition.list,
      learning.systems.list,
      instructional.design.list,
      field.experiences.list,
      ela.list,
      educational.change.list,
      attitudes.beliefs.list,
      scale.social.tech.list,
      scale.subjects.list,
      scale.tech.in.field.list,
      scale.phenomenology.list,
      scale.methods.list,
      scale.pipeline.list,
      scale.justice.list,
      scale.curriculum.list,
      social.scale.list
    )
  )
list.of.lists <- as.data.frame(list.of.lists)
write_csv(list.of.lists, "output/all_descriptors.csv")

# Bring the data back into a dataframe and return the key manuscriptIDs to row names
analysis.scales <- as.data.frame(analysis.scales)
rownames(analysis.scales) <- analysis.scales$manuscriptID
analysis.scales$manuscriptID <- NULL

# Scale the values so they are continuous values and normalized across the dataframe
analysis.scales <- scale(analysis.scales)

# Rename columns to easier-to-read descriptions
colnames(analysis.scales) <-
  c(
    "instructional_design",
    "specific_technologies",
    "attitudes_beliefs",
    "simulations_games",
    "teacher_pipeline",
    "methods",
    "phenomonology",
    "subject_areas",
    "socializing_technology_use"
  )

analysis.scales <- as.data.frame(analysis.scales)
analysis.scales$manuscriptID <- rownames(analysis.scales)
rownames(analysis.scales) <- NULL

#write_csv(as.data.frame(analysis.scales), "data/analysis-full.csv")



############################################################
#  Compute hierarchical clustering on principal components #
############################################################

analysis.scales <-
  read_csv("data/analysis-full.csv",
           col_names = TRUE,
           show_col_types = FALSE)
analysis.scales <- as.data.frame(analysis.scales)
rownames(analysis.scales) <- analysis.scales$manuscriptID
analysis.scales$manuscriptID <- NULL

# Re-run PCA on newly compiled and scaled dataframe
analysis.pca <- PCA(analysis.scales, graph = FALSE)

# Calculate HCPC
res.hcpc <- HCPC(analysis.pca, consol = FALSE, graph = FALSE)

analysis.cluster.frame <<- as.data.frame(res.hcpc$data.clust$clust)
analysis.cluster.frame$manuscriptID <- rownames(res.hcpc$data.clust)
analysis.cluster.frame <- analysis.cluster.frame %>%
  rename("cluster" = "res.hcpc$data.clust$clust") %>%
  select(manuscriptID, cluster)

analysis.cluster.frame$color <-
  iu.colors.fade[analysis.cluster.frame$cluster]

#################
# Create Graphs #
#################

# Create and output dendrogram of clusters by manuscriptID
fviz_dend(
  res.hcpc,
  cex = 0.7,
  # Label size
  palette = iu.colors,
  # Color palette see ?ggpubr::ggpar
  rect = FALSE,
  rect_fill = FALSE,
  # Add rectangle around groups
  ggtheme = theme_minimal(),
  labels_track_height = 0.8
) +     # Augment the room for labels
  scale_x_continuous(breaks = NULL)
ggsave(
  "output/cluster_dend.pdf",
  width = 22,
  height = 6,
  units = "in",
  dpi = 300
)

# Create and output a cluster map using calculated X-Y coordinates
fviz_cluster(
  res.hcpc,
  repel = TRUE,
  # Avoid label overlapping
  show.clust.cent = TRUE,
  # Show cluster centers
  palette = iu.colors,
  # Color palette see ?ggpubr::ggpar
  ggtheme = theme_minimal(),
  axes = c(1, 2),
  ylab = "Attitudes and Beliefs (-) - Behaviors and Practices (+)",
  xlab = "Teacher Pipeline (-) - Social and Individual Experiences (+)",
  main = ""
)
ggsave(
  "output/cluster_map1-2.pdf",
  width = 11,
  height = 8.5,
  units = "in",
  dpi = 300,
  bg = "#FFFFFF"
)

# Create and output a cluster map using calculated X-Y coordinates
fviz_cluster(
  res.hcpc,
  repel = TRUE,
  # Avoid label overlapping
  show.clust.cent = TRUE,
  # Show cluster centers
  palette = iu.colors,
  # Color palette see ?ggpubr::ggpar
  ggtheme = theme_minimal(),
  axes = c(1, 3),
  main = "Cluster Map of Calculated Clusters"
)
ggsave(
  "output/cluster_map1-3.pdf",
  width = 11,
  height = 8.5,
  units = "in",
  dpi = 300,
  bg = "#FFFFFF"
)

# Create and output a cluster map using calculated X-Y coordinates
fviz_cluster(
  res.hcpc,
  repel = TRUE,
  # Avoid label overlapping
  show.clust.cent = TRUE,
  # Show cluster centers
  palette = iu.colors,
  # Color palette see ?ggpubr::ggpar
  ggtheme = theme_minimal(),
  axes = c(2, 3),
  main = "Cluster Map of Calculated Clusters"
)
ggsave(
  "output/cluster_map2-3.pdf",
  width = 11,
  height = 8.5,
  units = "in",
  dpi = 300,
  bg = "#FFFFFF"
)


# Create quadrant chart
quadrant.values <- c(16, 24, 26, 21)
quadrant.plot <- four_quadrant(quadrant.values)
quadrant.plot <- quadrant.plot +
  xlim(-10, 10) + ylim(-5.4, 5.4) +
  annotate(
    geom = "text",
    x = 5,
    y = -2,
    label = "Studies that explore\nsocial and individual experiences\nthrough attitudes and beliefs",
    hjust = 0,
    vjust = 1,
    size = 3.5,
    color = "#243142"
  ) +
  annotate(
    geom = "text",
    x = 4.5,
    y = 2.5,
    label = "Studies that explore\nsocial and individual experiences\nthrough behaviors and practices",
    hjust = 0,
    vjust = 1,
    size = 3.5,
    color = "#243142"
  ) +
  annotate(
    geom = "text",
    x = -5.5,
    y = -2.5,
    label = "Studies that explore\nthe teacher pipeline\nthrough attitudes and beliefs",
    hjust = 1,
    vjust = 1,
    size = 3.5,
    color = "#243142"
  ) +
  annotate(
    geom = "text",
    x = -5.25,
    y = 3,
    label = "Studies that explore\nthe teacher pipeline\nthrough behaviors and practices",
    hjust = 1,
    vjust = 1,
    size = 3.5,
    color = "#243142"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  )
quadrant.plot
ggsave(
  "output/quadrant-chart.pdf",
  plot = quadrant.plot,
  width = 11,
  height = 4.25,
  units = "in",
  dpi = 300,
  bg = "#FFFFFF"
)

# Create and output scree plot to identify explanatory dimensions
fviz_eig(
  analysis.pca,
  addlabels = TRUE,
  ylim = c(0, 20),
  ggtheme = theme_minimal()
) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())
ggsave(
  "output/scree_plot.pdf",
  width = 11,
  height = 8.5,
  units = "in",
  dpi = 300,
  bg = "#FFFFFF"
)

# Pull out variables for analysis
var <- get_pca_var(analysis.pca)

# Create and output a correlation plot for factors and dimensions. Of particular
# interest are dimensions 1 and 2, and potentially 3.
ggcorrplot(t(var$cos2),
           method = "circle",
           colors = c("#990000", "#59264D", "#FFAA00")) +
  labs(x = "Dimensions", y = "Factors", title = "Dimension Correlations")
ggsave(
  "output/corrplot.pdf",
  width = 8.5,
  height = 11,
  units = "in",
  dpi = 300
)

# Scree plots of each dimension as another representation of the correlation plot
# Dimension 1
fviz_contrib(analysis.pca,
             choice = "var",
             axes = 1,
             top = 10)
ggsave(
  "output/dim1-contributions.pdf",
  width = 11,
  height = 8.5,
  units = "in",
  dpi = 300
)
# Dimension 2
fviz_contrib(analysis.pca,
             choice = "var",
             axes = 2,
             top = 10)
ggsave(
  "output/dim2-contributions.pdf",
  width = 11,
  height = 8.5,
  units = "in",
  dpi = 300
)
# Dimension 3
fviz_contrib(analysis.pca,
             choice = "var",
             axes = 3,
             top = 10)
ggsave(
  "output/dim3-contributions.pdf",
  width = 11,
  height = 8.5,
  units = "in",
  dpi = 300
)


# Dim 1: Teacher Pipeline (-) - Social and Individual Experiences (+)
# Dim 2: Attitudes and Beliefs (-) - Behaviors and Practices (+)

# PCA factor map color by cos2 values: quality on the factor map
fviz_pca_var(
  analysis.pca,
  col.var = "cos2",
  gradient.cols = c("#59264D", "#FFAA00", "#990000"),
  repel = TRUE,
  # Avoid text overlapping
)
ggsave(
  "output/pca.pdf",
  width = 11,
  height = 11,
  units = "in",
  dpi = 300
)



# Identify the "paragons" from the HCPC analysis, that is, those papers that are
# closest to the center of each cluster, indicating they are good representatives.
res.hcpc$desc.ind$para

# Pull out contributing factor details from each cluster and turn it into a dataframe

# Cluster 1
cluster.details.1 <-
  frame.factor.details(as.data.frame(res.hcpc$desc.var$quanti$`1`), 1)
# Cluster 1: Learning to Teach with Games and Simulations

# Cluster 2
cluster.details.2 <-
  frame.factor.details(as.data.frame(res.hcpc$desc.var$quanti$`2`), 2)
# Cluster 2: Building a TPACK-Informed Teacher Workforce

# Cluster 3
cluster.details.3 <-
  frame.factor.details(as.data.frame(res.hcpc$desc.var$quanti$`3`), 3)
# Cluster 3: Understanding What Teachers Do with Technology

# Cluster 4
cluster.details.4 <-
  frame.factor.details(as.data.frame(res.hcpc$desc.var$quanti$`4`), 4)
# Cluster 4: Shaping Beliefs and Attitudes around Technology Use through Socialization

# Cluster 5
cluster.details.5 <-
  frame.factor.details(as.data.frame(res.hcpc$desc.var$quanti$`5`), 5)
# Understanding Impressions Around Technology Use in the Field

# Bind all the clusters together into one dataframe for easier use and display
cluster.details.frame <-
  bind_rows(
    cluster.details.1,
    cluster.details.2,
    cluster.details.3,
    cluster.details.4,
    cluster.details.5
  )

write_csv(cluster.details.frame, "output/cluster_details.csv")

factor.flow.frame <-
  read_csv(
    "data/factor_flow-no_exclude.csv",
    col_names = TRUE,
    show_col_types = FALSE
  )
ggplot(data = factor.flow.frame,
       aes(
         axis1 = start,
         axis2 = descriptor_reduction,
         axis3 = final
       )) +
  scale_x_discrete(limits = c("Initial", "PCA Reduction", #"Scale Development",
                              "Final"),
                   expand = c(.2, .05)) +
  xlab("Process") +
  geom_alluvium(aes(fill = final)) +
  geom_stratum(alpha = 0, color = "#EDEBEB36") +
  ggfittext::geom_fit_text(
    stat = "stratum",
    aes(label = after_stat(stratum)),
    width = 1 / 4,
    min.size = 2.5,
    color = "#191919"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank()
  ) +
  scale_fill_manual(values = iu.10.colors) #+
ggsave(
  "output/process-alluvial.pdf",
  width = 22,
  height = 11,
  units = "in",
  dpi = 300
)


compare.count <-
  read_csv("data/te-tet-differences.csv",
           col_names = TRUE,
           show_col_types = FALSE)
compare.count <- subset(compare.count,!(terms %in% "baseline"))
compare.count$te2 <-
  round((compare.count$te + 1) / sum(compare.count$te + 1), digits = 3)
compare.count$tet2 <-
  round((compare.count$tet + 1) / sum(compare.count$tet + 1), digits = 3)
compare.count$tet2[compare.count$tet == 0] <- 0
compare.count$difference <-
  round(log2(compare.count$te2 / compare.count$tet2), digits = 3)
compare.count$difference[compare.count$difference > 3] <- 3
compare.count$difference[compare.count$difference < -3] <- -3
compare.count[compare.count$difference > 0, "status"] <-
  "Teacher Education"
compare.count[compare.count$difference < 0, "status"] <-
  "Teacher Education with Technology"
library(ggpubr)
compare.graph2 <-
  ggdotchart(
    compare.count,
    x = "terms",
    y = "difference",
    color = "status",
    sorting = "descending",
    add = "segments",
    rotate = TRUE,
    dot.size = 6,
    palette = iu.colors,
    # label = abs(round(compare.count$difference, 1)),
    #font.label = list(size = 8, vjust = 0.5, color = "white"),
    ggtheme = theme_pubr()
  ) +
  geom_hline(yintercept = 0,
             linetype = 3,
             color = "lightgray") +
  scale_y_continuous(
    breaks = -3:3,
    limits = c(-3, 3),
    labels = c("8X+", "4X", "2X", "Same", "2X", "4X", "8X+")
  ) + xlab("Descriptors") +
  labs(color = "") +
  ylab("Relative difference (log2)") +
  theme(text = element_text(size = 10))
compare.graph2
ggsave(
  "output/te-tet-compare.pdf",
  width = 8.5,
  height = 11,
  units = "in",
  dpi = 300
)

#################################
# Network analysis of citations #
#################################

analysis.cites <- read_csv(
  "data/mtf-citations-edgelist.csv",
  col_names = TRUE,
  show_col_types = FALSE
)
analysis.cites <-
  analysis.cites[analysis.cites$cite %in% names(which(table(analysis.cites$cite) > 1)),]

# OK, move through different clusters...
# Full Graph
cites.all.centrality <- calculate.centrality(analysis.cites)
cites.6.keyactors <- calculate.keyactors(cites.all.centrality, 6)
network.plot(analysis.cites, cites.6.keyactors, 0, 4, 50, 22, 17)
write_csv(cites.6.keyactors, file = "output/keyactors-all.csv", append = FALSE)


# Cluster 1
cites.1.frame <- isolate.cluster(1, analysis.cites)
cites.1.centrality <- calculate.centrality(cites.1.frame)
cites.1.keyactors <- calculate.keyactors(cites.1.centrality, 1)
network.plot(cites.1.frame, cites.1.keyactors, 1, 4, 20, 11, 8.5)

# Cluster 2
cites.2.frame <- isolate.cluster(2, analysis.cites)
cites.2.centrality <- calculate.centrality(cites.2.frame)
cites.2.keyactors <- calculate.keyactors(cites.2.centrality, 2)
network.plot(cites.2.frame, cites.2.keyactors, 2, 4, 40, 22, 17)

# Cluster 3
cites.3.frame <- isolate.cluster(3, analysis.cites)
cites.3.centrality <- calculate.centrality(cites.3.frame)
cites.3.keyactors <- calculate.keyactors(cites.3.centrality, 3)
network.plot(cites.3.frame, cites.3.keyactors, 3, 2, 40, 11, 8.5)

# Cluster 4
cites.4.frame <- isolate.cluster(4, analysis.cites)
cites.4.centrality <- calculate.centrality(cites.4.frame)
cites.4.keyactors <- calculate.keyactors(cites.4.centrality, 4)
network.plot(cites.4.frame, cites.4.keyactors, 4, 2, 40, 11, 8.5)

# Cluster 5
cites.5.frame <- isolate.cluster(5, analysis.cites)
cites.5.centrality <- calculate.centrality(cites.5.frame)
cites.5.keyactors <- calculate.keyactors(cites.5.centrality, 5)
network.plot(cites.5.frame, cites.5.keyactors, 5, 1, 20, 11, 8.5)

cites.keyactors <- bind_rows(
  cites.1.keyactors,
  cites.2.keyactors,
  cites.3.keyactors,
  cites.4.keyactors,
  cites.5.keyactors
) %>%
  arrange(cluster,
          desc(keystatus),
          desc(interaction_score),
          manuscriptID)
write_csv(cites.keyactors, file = "output/keyactors.csv", append = FALSE)

