import groovy.json.JsonSlurper
import java.text.SimpleDateFormat
// Used for checking versions of nextclade
def compareVersions(String version1, String version2)
{
    if (version2 == null){
      return 0;
    }
    String[] string1Vals = version1.split("\\.");
    String[] string2Vals = version2.split("\\.");

    int length = Math.max(string1Vals.length, string2Vals.length);

    for (int i = 0; i < length; i++)
    {
        Integer v1 = (i < string1Vals.length)?Integer.parseInt(string1Vals[i]):0;
        Integer v2 = (i < string2Vals.length)?Integer.parseInt(string2Vals[i]):0;

        //Making sure Version1 bigger than version2
        if (v1 > v2)
        {
            return 1;
        }
        //Making sure Version1 smaller than version2
        else if(v1 < v2)
        {
            return -1;
        }
    }

    //Both are equal
    return 0;
}

def checkVersions(nextclade_version,dir)
{
  jsonSlurper = new JsonSlurper()
  tagData = jsonSlurper.parseText(file(dir.resolve('./files/tag.json')).text)

  if(
      compareVersions(params.nextclade_version,tagData.compatibility.nextcladeCli.min) != -1 &&
      compareVersions(params.nextclade_version,tagData.compatibility.nextcladeCli.max) != 1
  ){
     // ok this directory is >= the min and <= the max version supported so add to our list of compatible data
     return true;
  }else{
    return false;
  }


}


def nextcladeVersionChecker(nextclade_data_tag){
  if (nextclade_data_tag == null) {

    tagged_dirs = file(projectDir.resolve("./data/nextclade/datasets/sars-cov-2/references/MN908947/versions/*"), type: 'dir', maxdepth: 1)

    // discard those not compatible with current selected version
    def compatible_tagged_dirs = []
    for (dir in tagged_dirs){
      if(checkVersions(params.nextclade_version,dir)){
        compatible_tagged_dirs.add(dir)
      }
    }

    // now we have all the compatible data tags get the most recent
    def tag_list = []
    date_parse = new SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss'Z'");
    compatible_tagged_dirs.each { val -> tag_list << date_parse.parse(val.getBaseName()) }
    nextclade_data_tag = Collections.max(tag_list).format("yyyy-MM-dd'T'HH:mm:ss'Z'")

  } else{

    tagDir = file(projectDir.resolve("./data/nextclade/datasets/sars-cov-2/references/MN908947/versions/").resolve(nextclade_data_tag))

    if(checkVersions(params.nextclade_version,tagDir)){
      println("Nextclade $params.nextclade_version is compatible with $nextclade_data_tag")
    } else{

      //dataset is incompatible
      println("Nextclade $params.nextclade_version is incompatible with $nextclade_data_tag")
      exit 1

    }

  }

  nextclade = projectDir.resolve('./data/nextclade/datasets/sars-cov-2/references/MN908947/versions/').resolve(nextclade_data_tag).resolve('./files')
  nextclade_dataset = file(nextclade, type: 'dir', checkIfExists:true)
  return [nextclade_data_tag,nextclade_dataset]
}
