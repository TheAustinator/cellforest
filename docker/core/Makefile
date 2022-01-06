include make_env

.PHONY: build rebuild shell connect stop rm

build: Dockerfile
	docker build $(BUILD_ARGS) -t $(NS)/$(IMAGE_NAME):$(VERSION) -f Dockerfile .

rebuild: Dockerfile
	docker build --no-cache -t $(NS)/$(IMAGE_NAME):$(VERSION) -f Dockerfile .

shell:
	docker run --rm --name $(CONTAINER_NAME)-$(CONTAINER_INSTANCE) -it $(PORTS) $(VOLUMES) $(RUN_ARGS) $(NS)/$(IMAGE_NAME):$(VERSION)

connect:
	docker exec -it $(CONTAINER_NAME)-$(CONTAINER_INSTANCE) zsh

stop:
	docker stop $(CONTAINER_NAME)-$(CONTAINER_INSTANCE)

rm:
	docker rm $(CONTAINER_NAME)-$(CONTAINER_INSTANCE)

default: build
