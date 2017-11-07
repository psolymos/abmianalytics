library(sendmailR)
from <- sprintf("<WebSI_Team@\\%s>", Sys.info()[4])
to <- "<psolymos@gmail.com>"
subject <- "Your intactness results are ready"
## change iris to the zip file
body <- list("Hi,\n\nYour intactness results are ready, see attachment.\n\nWith regards,\n\nthe WebSI Team", mime_part(iris))
sendmail(from, to, subject, body,
         control=list(smtpServer="ASPMX.L.GOOGLE.COM"))
