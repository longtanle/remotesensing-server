from bson.objectid import ObjectId

from main.db import MongoDB


class EarthExplorerService:
    """ doc string for EarthExplorerService """

    def __init__(self):
        super(EarthExplorerService, self).__init__()
        self.collection = "earthexplorers"
        self.mongo = MongoDB()

    def add(self, earthexplorer_obj):
        earthexplorer = self.mongo.find(self.collection, {"book_name": earthexplorer_obj["book_name"]})
        if not earthexplorer:
            return (
                self.mongo.save(self.collection, earthexplorer_obj),
                "Successfully created.",
                200,
            )
        else:
            return ("ok", "Book already added to the library.", 400)

    def earthexplorers_list(self):
        return self.mongo.find(self.collection)

    def delete_earthexplorer(self, earthexplorer_id):
        return self.mongo.delete(self.collection, earthexplorer_id)

    def update_earthexplorer(self, earthexplorer_id, earthexplorer_obj):
        condition = {"$set": earthexplorer_obj}
        res, update_count = self.mongo.update(self.collection, earthexplorer_id, condition)

        if res:
            return ("success", res, "ok", 200)
        return ("error", "", "Something went wrong.", 400)

    def get_earthexplorer(self, earthexplorer_id):
        condition = {"_id": ObjectId(earthexplorer_id)}
        return self.mongo.find(self.collection, condition)
