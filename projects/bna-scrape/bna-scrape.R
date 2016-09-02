library(XML)
library(RCurl)
library(mefa4)

scrapeBNA <- function(species) {
    url <- paste0("https://www.allaboutbirds.org/guide/", species, "/lifehistory/index.html")
    page <- getURL(url, header=FALSE, verbose=TRUE)
    pagetree <- htmlTreeParse(page, error=function(...){}, useInternalNodes = TRUE)
    parseTree <- function(pagetree, id) {
        x <- getNodeSet(pagetree, paste0("//*/div[@id='", id, "']"))
        x <- paste(capture.output(print(x[[1]])), collapse = "")
        x <- strsplit(x, "p>")[[1]][2]
        substr(x, 1,nchar(x)-2)
    }
    list(habitat = parseTree(pagetree, "life_habitat"),
        food = parseTree(pagetree, "life_food"),
        nesting = parseTree(pagetree, "life_nesting"),
        behavior = parseTree(pagetree, "life_behavior"))
}

x1 <- getNodeSet(pagetree, paste0("//*/div[@id='life_habitat']"))
x2 <- getNodeSet(pagetree, paste0("//*/div[@id='life_food']"))
x3 <- getNodeSet(pagetree, paste0("//*/div[@id='life_nesting']"))
x4 <- getNodeSet(pagetree, paste0("//*/div[@id='life_behavior']"))

scrapeBNA("Connecticut_Warbler")
scrapeBNA("Acadian_Flycatcher")

all_spp <- c("Abert's Towhee", "Acadian Flycatcher", "Acorn Woodpecker",
    "African Collared-Dove", "Alder Flycatcher", "Allen's Hummingbird",
    "Altamira Oriole", "American Avocet", "American Bittern", "American Black Duck",
    "American Coot", "American Crow", "American Dipper", "American Golden-Plover",
    "American Goldfinch", "American Kestrel", "American Oystercatcher",
    "American Pipit", "American Redstart", "American Robin", "American Three-toed Woodpecker",
    "American Tree Sparrow", "American White Pelican", "American Wigeon",
    "American Woodcock", "Ancient Murrelet", "Anhinga", "Anna's Hummingbird",
    "Arctic Tern", "Arctic Warbler", "Arizona Woodpecker", "Ash-throated Flycatcher",
    "Atlantic Puffin", "Audubon's Oriole", "Baird's Sandpiper", "Baird's Sparrow",
    "Bald Eagle", "Baltimore Oriole", "Band-tailed Pigeon", "Bank Swallow",
    "Barn Owl", "Barn Swallow", "Barred Owl", "Barrow's Goldeneye",
    "Bay-breasted Warbler", "Bell's Sparrow", "Bell's Vireo", "Belted Kingfisher",
    "Bendire's Thrasher", "Bewick's Wren", "Bicknell's Thrush", "Black-and-white Warbler",
    "Black-backed Woodpecker", "Black-bellied Plover", "Black-bellied Whistling-Duck",
    "Black-billed Cuckoo", "Black-billed Magpie", "Black-capped Chickadee",
    "Black-capped Vireo", "Black-chinned Hummingbird", "Black-crested Titmouse",
    "Black-crowned Night-Heron", "Black-footed Albatross", "Black-headed Grosbeak",
    "Black-headed Gull", "Black-legged Kittiwake", "Black-necked Stilt",
    "Black-tailed Gnatcatcher", "Black-throated Blue Warbler", "Black-throated Gray Warbler",
    "Black-throated Green Warbler", "Black-throated Sparrow", "Black-whiskered Vireo",
    "Black Guillemot", "Black Oystercatcher", "Black Phoebe", "Black Rail",
    "Black Rosy-Finch", "Black Scoter", "Black Skimmer", "Black Tern",
    "Black Turnstone", "Black Vulture", "Blackburnian Warbler", "Blackpoll Warbler",
    "Blue-footed Booby", "Blue-gray Gnatcatcher", "Blue-headed Vireo",
    "Blue-throated Hummingbird", "Blue-winged Teal", "Blue-winged Warbler",
    "Blue Grosbeak", "Blue Jay", "Bluethroat", "Boat-tailed Grackle",
    "Bobolink", "Bohemian Waxwing", "Bonaparte's Gull", "Boreal Chickadee",
    "Boreal Owl", "Botteri's Sparrow", "Brandt's Cormorant", "Brant",
    "Brewer's Blackbird", "Brewer's Sparrow", "Bridled Titmouse",
    "Broad-billed Hummingbird", "Broad-tailed Hummingbird", "Broad-winged Hawk",
    "Bronzed Cowbird", "Brown-capped Rosy-Finch", "Brown-headed Cowbird",
    "Brown-headed Nuthatch", "Brown Booby", "Brown Creeper", "Brown Pelican",
    "Brown Thrasher", "Buff-bellied Hummingbird", "Buff-breasted Flycatcher",
    "Buff-breasted Sandpiper", "Bufflehead", "Bullock's Oriole",
    "Burrowing Owl", "Bushtit", "Cackling Goose", "Cactus Wren",
    "California Condor", "California Gull", "California Quail", "California Scrub-Jay",
    "California Thrasher", "California Towhee", "Calliope Hummingbird",
    "Canada Goose", "Canada Warbler", "Canvasback", "Canyon Towhee",
    "Canyon Wren", "Cape May Warbler", "Carolina Chickadee", "Carolina Wren",
    "Caspian Tern", "Cassin's Auklet", "Cassin's Finch", "Cassin's Kingbird",
    "Cassin's Sparrow", "Cassin's Vireo", "Cattle Egret", "Cave Swallow",
    "Cedar Waxwing", "Cerulean Warbler", "Chestnut-backed Chickadee",
    "Chestnut-collared Longspur", "Chestnut-sided Warbler", "Chihuahuan Raven",
    "Chimney Swift", "Chipping Sparrow", "Chuck-will's-widow", "Chukar",
    "Cinnamon Teal", "Clapper Rail", "Clark's Grebe", "Clark's Nutcracker",
    "Clay-colored Sparrow", "Cliff Swallow", "Common Eider", "Common Gallinule",
    "Common Goldeneye", "Common Grackle", "Common Ground-Dove", "Common Loon",
    "Common Merganser", "Common Murre", "Common Nighthawk", "Common Pauraque",
    "Common Poorwill", "Common Raven", "Common Redpoll", "Common Tern",
    "Common Yellowthroat", "Connecticut Warbler", "Cooper's Hawk",
    "Cordilleran Flycatcher", "Costa's Hummingbird", "Couch's Kingbird",
    "Crested Caracara", "Crissal Thrasher", "Curve-billed Thrasher",
    "Dark-eyed Junco", "Dickcissel", "Double-crested Cormorant",
    "Dovekie", "Downy Woodpecker", "Dunlin", "Dusky Flycatcher",
    "Dusky Grouse", "Eared Grebe", "Eastern Bluebird", "Eastern Kingbird",
    "Eastern Meadowlark", "Eastern Phoebe", "Eastern Screech-Owl",
    "Eastern Towhee", "Eastern Whip-poor-will", "Eastern Wood-Pewee",
    "Elegant Tern", "Elegant Trogon", "Emperor Goose", "Eurasian Collared-Dove",
    "Eurasian Tree Sparrow", "Eurasian Wigeon", "European Starling",
    "Evening Grosbeak", "Ferruginous Hawk", "Field Sparrow", "Fish Crow",
    "Flammulated Owl", "Florida Scrub-Jay", "Fork-tailed Storm-Petrel",
    "Forster's Tern", "Fox Sparrow", "Franklin's Gull", "Fulvous Whistling-Duck",
    "Gadwall", "Gambel's Quail", "Gila Woodpecker", "Gilded Flicker",
    "Glaucous-winged Gull", "Glaucous Gull", "Glossy Ibis", "Golden-cheeked Warbler",
    "Golden-crowned Kinglet", "Golden-crowned Sparrow", "Golden-fronted Woodpecker",
    "Golden-winged Warbler", "Golden Eagle", "Grace's Warbler", "Grasshopper Sparrow",
    "Gray-cheeked Thrush", "Gray-crowned Rosy-Finch", "Gray Catbird",
    "Gray Flycatcher", "Gray Hawk", "Gray Jay", "Gray Partridge",
    "Gray Vireo", "Great-tailed Grackle", "Great Black-backed Gull",
    "Great Blue Heron", "Great Cormorant", "Great Crested Flycatcher",
    "Great Egret", "Great Gray Owl", "Great Horned Owl", "Great Kiskadee",
    "Greater Pewee", "Greater Prairie-Chicken", "Greater Roadrunner",
    "Greater Sage-Grouse", "Greater Scaup", "Greater White-fronted Goose",
    "Greater Yellowlegs", "Green-tailed Towhee", "Green-winged Teal",
    "Green Heron", "Green Jay", "Groove-billed Ani", "Gull-billed Tern",
    "Gunnison Sage-Grouse", "Gyrfalcon", "Hairy Woodpecker", "Hammond's Flycatcher",
    "Harlequin Duck", "Harris's Hawk", "Harris's Sparrow", "Heermann's Gull",
    "Henslow's Sparrow", "Hepatic Tanager", "Hermit Thrush", "Hermit Warbler",
    "Herring Gull", "Hoary Redpoll", "Hooded Merganser", "Hooded Oriole",
    "Hooded Warbler", "Horned Grebe", "Horned Lark", "Horned Puffin",
    "House Finch", "House Sparrow", "House Wren", "Hudsonian Godwit",
    "Hutton's Vireo", "Iceland Gull", "Inca Dove", "Indigo Bunting",
    "Ivory-billed Woodpecker", "Ivory Gull", "Juniper Titmouse",
    "Kentucky Warbler", "Killdeer", "King Eider", "King Rail", "Kirtland's Warbler",
    "Ladder-backed Woodpecker", "Lapland Longspur", "Lark Bunting",
    "Lark Sparrow", "Laughing Gull", "Lawrence's Goldfinch", "Laysan Albatross",
    "Lazuli Bunting", "Le Conte's Sparrow", "Le Conte's Thrasher",
    "Least Bittern", "Least Flycatcher", "Least Grebe", "Least Sandpiper",
    "Least Tern", "Lesser Black-backed Gull", "Lesser Goldfinch",
    "Lesser Prairie-Chicken", "Lesser Scaup", "Lesser Yellowlegs",
    "Lewis's Woodpecker", "Limpkin", "Lincoln's Sparrow", "Little Blue Heron",
    "Little Gull", "Loggerhead Shrike", "Long-billed Curlew", "Long-billed Dowitcher",
    "Long-billed Thrasher", "Long-eared Owl", "Long-tailed Duck",
    "Louisiana Waterthrush", "Lucifer Hummingbird", "Lucy's Warbler",
    "MacGillivray's Warbler", "Magnificent Frigatebird", "Magnificent Hummingbird",
    "Magnolia Warbler", "Mallard", "Mangrove Cuckoo", "Marbled Godwit",
    "Marbled Murrelet", "Marsh Wren", "Masked Booby", "McCown's Longspur",
    "Merlin", "Mew Gull", "Mexican Jay", "Mississippi Kite", "Monk Parakeet",
    "Montezuma Quail", "Mottled Duck", "Mountain Bluebird", "Mountain Chickadee",
    "Mountain Plover", "Mountain Quail", "Mourning Dove", "Mourning Warbler",
    "Muscovy Duck", "Mute Swan", "Nashville Warbler", "Nelson's Sparrow",
    "Neotropic Cormorant", "Northern Bobwhite", "Northern Cardinal",
    "Northern Flicker", "Northern Fulmar", "Northern Gannet", "Northern Goshawk",
    "Northern Harrier", "Northern Hawk Owl", "Northern Mockingbird",
    "Northern Parula", "Northern Pintail", "Northern Pygmy-Owl",
    "Northern Rough-winged Swallow", "Northern Saw-whet Owl", "Northern Shoveler",
    "Northern Shrike", "Northern Waterthrush", "Northwestern Crow",
    "Nuttall's Woodpecker", "Oak Titmouse", "Olive-sided Flycatcher",
    "Orange-crowned Warbler", "Orchard Oriole", "Osprey", "Ovenbird",
    "Pacific-slope Flycatcher", "Pacific Golden-Plover", "Pacific Loon",
    "Pacific Wren", "Painted Bunting", "Painted Redstart", "Palm Warbler",
    "Parakeet Auklet", "Pectoral Sandpiper", "Pelagic Cormorant",
    "Peregrine Falcon", "Phainopepla", "Philadelphia Vireo", "Pied-billed Grebe",
    "Pigeon Guillemot", "Pileated Woodpecker", "Pine Grosbeak", "Pine Siskin",
    "Pine Warbler", "Pinyon Jay", "Piping Plover", "Plain Chachalaca",
    "Plumbeous Vireo", "Prairie Falcon", "Prairie Warbler", "Prothonotary Warbler",
    "Purple Finch", "Purple Gallinule", "Purple Martin", "Purple Sandpiper",
    "Pygmy Nuthatch", "Pyrrhuloxia", "Razorbill", "Red-bellied Woodpecker",
    "Red-breasted Merganser", "Red-breasted Nuthatch", "Red-breasted Sapsucker",
    "Red-cockaded Woodpecker", "Red-eyed Vireo", "Red-faced Warbler",
    "Red-footed Booby", "Red-headed Woodpecker", "Red-naped Sapsucker",
    "Red-necked Grebe", "Red-shouldered Hawk", "Red-tailed Hawk",
    "Red-throated Loon", "Red-winged Blackbird", "Red Crossbill",
    "Red Knot", "Reddish Egret", "Redhead", "Rhinoceros Auklet",
    "Ridgway's Rail", "Ring-billed Gull", "Ring-necked Duck", "Ring-necked Pheasant",
    "Rock Pigeon", "Rock Ptarmigan", "Rock Sandpiper", "Rock Wren",
    "Rose-breasted Grosbeak", "Roseate Spoonbill", "Roseate Tern",
    "Ross's Goose", "Ross's Gull", "Rough-legged Hawk", "Royal Tern",
    "Ruby-crowned Kinglet", "Ruby-throated Hummingbird", "Ruddy Duck",
    "Ruddy Turnstone", "Ruffed Grouse", "Rufous-crowned Sparrow",
    "Rufous-winged Sparrow", "Rufous Hummingbird", "Rusty Blackbird",
    "Sabine's Gull", "Sage Thrasher", "Sagebrush Sparrow", "Saltmarsh Sparrow",
    "Sanderling", "Sandhill Crane", "Sandwich Tern", "Savannah Sparrow",
    "Say's Phoebe", "Scaled Quail", "Scarlet Tanager", "Scissor-tailed Flycatcher",
    "Scott's Oriole", "Seaside Sparrow", "Sedge Wren", "Semipalmated Plover",
    "Semipalmated Sandpiper", "Sharp-shinned Hawk", "Sharp-tailed Grouse",
    "Shiny Cowbird", "Short-billed Dowitcher", "Short-eared Owl",
    "Smith's Longspur", "Smooth-billed Ani", "Snail Kite", "Snow Bunting",
    "Snow Goose", "Snowy Egret", "Snowy Owl", "Snowy Plover", "Solitary Sandpiper",
    "Song Sparrow", "Sooty Grouse", "Sora", "Spot-breasted Oriole",
    "Spotted Owl", "Spotted Sandpiper", "Spotted Towhee", "Sprague's Pipit",
    "Spruce Grouse", "Steller's Eider", "Steller's Jay", "Stilt Sandpiper",
    "Summer Tanager", "Surf Scoter", "Surfbird", "Swainson's Hawk",
    "Swainson's Thrush", "Swainson's Warbler", "Swallow-tailed Kite",
    "Swamp Sparrow", "Tennessee Warbler", "Thayer's Gull", "Thick-billed Murre",
    "Townsend's Solitaire", "Townsend's Warbler", "Tree Swallow",
    "Tricolored Heron", "Tropical Kingbird", "Trumpeter Swan", "Tufted Puffin",
    "Tufted Titmouse", "Tundra Swan", "Turkey Vulture", "Upland Sandpiper",
    "Varied Bunting", "Varied Thrush", "Vaux's Swift", "Veery", "Verdin",
    "Vermilion Flycatcher", "Vesper Sparrow", "Violet-green Swallow",
    "Virginia Rail", "Wandering Tattler", "Warbling Vireo", "Western Bluebird",
    "Western Grebe", "Western Gull", "Western Kingbird", "Western Meadowlark",
    "Western Sandpiper", "Western Screech-Owl", "Western Tanager",
    "Western Wood-Pewee", "Whimbrel", "White-breasted Nuthatch",
    "White-crowned Pigeon", "White-crowned Sparrow", "White-eyed Vireo",
    "White-faced Ibis", "White-headed Woodpecker", "White-rumped Sandpiper",
    "White-tailed Hawk", "White-tailed Kite", "White-tailed Ptarmigan",
    "White-throated Sparrow", "White-throated Swift", "White-tipped Dove",
    "White-winged Crossbill", "White-winged Dove", "White-winged Scoter",
    "White Ibis", "Whooping Crane", "Wild Turkey", "Willet", "Williamson's Sapsucker",
    "Willow Flycatcher", "Willow Ptarmigan", "Wilson's Phalarope",
    "Wilson's Plover", "Wilson's Snipe", "Wilson's Warbler", "Winter Wren",
    "Wood Duck", "Wood Stork", "Wood Thrush", "Woodhouse's Scrub-Jay",
    "Worm-eating Warbler", "Wrentit", "Yellow-bellied Flycatcher",
    "Yellow-bellied Sapsucker", "Yellow-billed Cuckoo", "Yellow-billed Magpie",
    "Yellow-breasted Chat", "Yellow-crowned Night-Heron", "Yellow-headed Blackbird",
    "Yellow-rumped Warbler", "Yellow-throated Vireo", "Yellow-throated Warbler",
    "Yellow Rail", "Yellow Warbler", "Zone-tailed Hawk")

all <- gsub("'", "", all_spp, fixed = TRUE)
all <- gsub(" ", "_", all, fixed = TRUE)

OK <- list()
err <- list()
for (i in all) {
    cat(i, " (", which(all == i), "/", length(all), ")", sep = "")
    flush.console()
    tmp <- try(scrapeBNA(i), silent = TRUE)
    if (!inherits(tmp, "try-error")) {
        OK[[i]] <- tmp
        cat(" ... OK\n")
    } else {
        err[[i]] <- tmp
        cat(" ... Error\n")
    }
}

tab <- data.frame(species=all_spp, t(sapply(OK, unlist)))
write.csv(tab, "bna-life-history-info.csv")
