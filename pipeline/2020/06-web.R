library(jsonlite)

x <- fromJSON("https://science.abmi.ca/results/reports/2020/images/index.json")$species
OUT <- "~/repos/species-web-dev/"

i <- 1
for (i in 1819:nrow(x)) {
    cat(i, "/", nrow(x), g, s, "\n")
    flush.console()
    g <- x$taxonid[i]
    s <- x$id[i]
    d <- try(fromJSON(paste0(
        'https://science.abmi.ca/results/reports/2020/images/', g, '/', s, '/index.json')))
    if (s == "Petunia.x.atkinsiana")
        d <- fromJSON('{"id":"Petunia.x.atkinsiana","scientific":"Petunia x atkinsiana","common":"Common Petunia","display":"Petunia x atkinsiana - Common Petunia","det":true,"useavailnorth":false,"useavailsouth":false,"coefnorth":false,"coefsouth":false,"map":false,"sectornorth":false,"sectorsouth":false,"idprev":"Petunia","idnext":"Phacelia.franklinii","taxonid":"vplants","taxonname":"Vascular plants"}')
    txt <- sprintf('<template>
      <spp :blob="blob"></spp>
    </template>
    <script>
    export default {
      head () {
        return {
          title: this.blob.display
        }
      },
      data () {
        return {
          blob: %s
        }
      }
    }
    </script>
    ', toJSON(d, auto_unbox=TRUE))

    if (!dir.exists(file.path(OUT, "pages", g, s)))
        dir.create(file.path(OUT, "pages", g, s))
    writeLines(txt, file.path(OUT, "pages", g, s, "index.vue"))
}




