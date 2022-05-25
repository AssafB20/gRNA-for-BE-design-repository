library(shiny)
library(reticulate)
source_python("BE_script1.py")

ui <- fluidPage(
  titlePanel("BE TargetFinder"),
  
  sidebarLayout(
    sidebarPanel(
      helpText("Find gRNA for your DNA sequence"),
      
      #textInput("text", label = h3("Enter your sequence here"), value = ""),
      textAreaInput("list", label = h3("Enter your sequence here"), value = ""),
      
      radioButtons("radio", label = h3("Choose type of Base Editor"),
                   choices = list("Adenine Base Editor" = 'a', "Cytidine Base Editor" = 'c'), 
                   selected = 'a'),
      fileInput("file", label = h3("Upload your sequences text file here"))
    ),
    
    mainPanel(
      img(src = "https://images.squarespace-cdn.com/content/v1/5f63d517b1e60e15744c391d/1601868540170-K05ETSW8PZNKCLL49GDP/ke17ZwdGBToddI8pDm48kO-4phTpo_CyN6HNWGRPW1oUqsxRUqqbr1mOJYKfIPR7LoDQ9mXPOjoJoqy81S2I8N_N4V1vUb5AoIIIbLZhVYy7Mythp_T-mtop-vrsUOmeInPi9iDjx9w8K4ZfjXt2do69x95e3VK-BlrAfDwAAjDFuf6SxFnKYaHSrOPlQbTDG6v6ULRah83RgHXAWD5lbQ/Main%25252Bimage%25252Bto%25252Bcrop.jpg", height = 300, width = 900),

      verbatimTextOutput("results")
      
    )
  )
)

server <- function(input, output) {
 input_seq <- reactive({
   #targetfinder_string(input$file, input$radio)
   targetfinder_string(input$list, input$radio)
 })
   
  output$results <- renderText({ 
    paste("pam, gRNA, editing window, number of bystanders", input_seq(), sep = "\n")
  })
 
 #input_file <- reactive({
  # targetfinder_file(input$file, input$radio)
 #}) 
 
  #output$results <- renderText({ 
   #paste("(pam, gRNA, editing window, number of bystanders):", input_file())
    
  #})
}

shinyApp(ui, server)