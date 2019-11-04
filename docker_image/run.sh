# pull latest tag from https://hub.docker.com/r/bgruening/galaxy-stable/tags
sudo docker run --privileged -i -t -p 8989:80 -p 8029:22 -p 9009:9002 \
-e "GALAXY_CONFIG_ADMIN_USERS=raphenar@mcmaster.ca,mcarthua@mcmaster.ca,jansem1@mcmaster.ca" \
-e "GALAXY_CONFIG_MASTER_API_KEY=w5CTPRGz6kqtW9BAAAfu3w42Xy" \
-e "GALAXY_CONFIG_BRAND=marcel_card_megares" -v /home/raphenar/galaxy_docker_images/marcel_card_megares/docker_image/galaxy_storage/:/export/ \
bgruening/galaxy-stable@sha256:827e648d336503eedf32faaba1ff8d61b12c11449f63dd6610492f5a1d1ebe75
