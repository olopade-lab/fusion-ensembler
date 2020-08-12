#! /bin/bash

function cleanup {
    DOCKER_CONTAINER=$(cat ${container_id_file})
    if [[ -n "${DOCKER_CONTAINER}" ]]
    then
        docker stop "${DOCKER_CONTAINER}"
        docker rm --force  "${DOCKER_CONTAINER}"
    fi
    rm -f "${container_id_file}"
}

trap cleanup EXIT

# create a file to which docker writes the container id
container_id_file=$(mktemp docker_container_.XXXXXX)
rm ${container_id_file}

docker run -d --cidfile ${container_id_file} --rm "${@}"

exit_status=$(docker wait $(cat ${container_id_file}))
exit ${exit_status}
