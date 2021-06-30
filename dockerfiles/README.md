# dockerfiles

cd ~/src/modPhred/dockerfiles/3.6.1
version=$(~/src/modPhred/run --version 2> /dev/null)
guppyversion=$(basename $PWD)
name=modphred-$guppyversion
echo $name:$version

docker build --pull -t lpryszcz/$name:$version .

docker tag lpryszcz/$name:$version lpryszcz/$name:latest

docker push lpryszcz/$name:$version && docker push lpryszcz/$name:latest

