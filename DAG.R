library(dagitty)
library(ggdag)
library(ggplot2)

dag <- dagitty('dag {
  BaboonPosition    [pos="0,2"]
  VervetPosition    [pos="0,4"]
  Terrain           [pos="2,1"]
  Vegetation        [pos="2,5"]
  Distance          [pos="3,3"]
  ViewshedClearance [pos="4.5,2"]
  LocalVisibility   [pos="4.5,4"]
  Detection         [pos="6,3"    latent]
  VervetMovement    [pos="8,3"    outcome]

  BaboonPosition   -> Distance
  VervetPosition   -> Distance
  BaboonPosition   -> ViewshedClearance
  VervetPosition   -> ViewshedClearance
  Terrain          -> ViewshedClearance
  Terrain          -> LocalVisibility
  Vegetation       -> ViewshedClearance
  Vegetation       -> LocalVisibility
  Distance         -> Detection
  ViewshedClearance -> Detection
  LocalVisibility  -> Detection
  Detection        -> VervetMovement
}')

dag_tidy <- tidy_dagitty(dag)

#  labels
dag_tidy$data$name_label <- dplyr::recode(dag_tidy$data$name,
                                          "BaboonPosition"    = "Baboon\nposition",
                                          "VervetPosition"    = "Vervet\nposition",
                                          "Terrain"           = "Terrain",
                                          "Vegetation"        = "Vegetation",
                                          "Distance"          = "Distance",
                                          "ViewshedClearance" = "Viewshed\nclearance",
                                          "LocalVisibility"   = "Local\nvisibility",
                                          "Detection"         = "Detection\n(latent)",
                                          "VervetMovement"    = "Vervet\nmovement"
)

# Colour by role
dag_tidy$data$node_type <- dplyr::case_when(
  dag_tidy$data$name %in% c("BaboonPosition", "VervetPosition") ~ "Position",
  dag_tidy$data$name %in% c("Terrain", "Vegetation") ~ "Landscape",
  dag_tidy$data$name %in% c("Distance", "ViewshedClearance", "LocalVisibility") ~ "Covariate",
  dag_tidy$data$name == "Detection" ~ "Latent",
  dag_tidy$data$name == "VervetMovement" ~ "Outcome"
)

p <- ggplot(dag_tidy, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_dag_edges(edge_width = 0.4, edge_colour = "grey40",
                 arrow_directed = grid::arrow(length = grid::unit(6, "pt"),
                                              type = "closed")) +
  geom_dag_point(aes(colour = node_type), size = 18) +
  geom_dag_text(aes(label = name_label), colour = "white", size = 2.5,
                lineheight = 0.9, fontface = "bold") +
  scale_colour_manual(
    values = c(Position = "#5C6BC0", Landscape = "#26A69A",
               Covariate = "#2E86AB", Latent = "#EF5350", Outcome = "#FF8F00"),
    name = NULL
  ) +
  theme_dag_blank() +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    legend.position = "bottom",
    legend.text = element_text(size = 10)
  )

p
ggsave("dag_clean.png", p, width = 11, height = 6, dpi = 300)